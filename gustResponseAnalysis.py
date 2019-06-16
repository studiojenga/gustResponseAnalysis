import csv, math, numpy
from datetime import datetime
from datetime import timedelta
import multiprocessing as mp

roh = 0.00012
alpha = 1.0 / 8.0
K1 = 8.0 #decay factor for x
K2 = 8.0 #decay factor for z
T = 600.0

class loc:
    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z

class dof:
    def __init__(self, x = 0.0, y = 0.0, z = 0.0, rx = 0.0, ry = 0.0, rz = 0.0):
        self.x = x
        self.y = y
        self.z = z
        self.rx = rx
        self.ry = ry
        self.rz = rz

class elm:
    def __init__(self, type = 'girder' , x = 0.0, y = 0.0, z = 0.0, l = 0.0, h = 0.0, b = 0.0, phai = 1.0, cD = 0.0, dcLda = 0.0):
        self.type = type
        self.loc = loc(x, y, z)
        self.l = l
        self.h = h
        self.b = b
        self.phai = phai
        self.cD = cD
        self.dcLda = dcLda
        self.mode = {}
        self.mass = dof()
        #used during calculation
        self.U = 0.0
        self.Iu = 0.0
        self.aa = {} #Aerodynamic admittance
        self.result = {}

class elm2:
    def __init__(self):
        self.corrF = 0.0
        self.deltaX = 0.0
        self.deltaZ = 0.0
        self.ps = {}
        self.z = 0.0
        self.U = 0.0
        self.Iu = 0.0

class mode:
    def __init__(self, omega = 0.0, dampD = 0.0, dampL = 0.0, dampM = 0.0, scale = 1.0):
        self.omega = omega
        self.damp = {'drag': dampD, 'lift': dampL, 'moment': dampM}
        self.scale = scale
        #used during calculation
        self.adamp = {}

class freq:
    def __init__(self, _freq, d_f = 0.0, psw = 0.0):
        self.freq = _freq
        self.d_f = d_f
        self.omega = 2.0 * math.pi * _freq
        #used during calculation
        self.psw = psw
        self.result = 0.0

def Uz(_z, _U10):
    return _U10 * (_z / 10.0) ** alpha

def Iuz(_z, _Iu10):
    return _Iu10 * (_z / 10.0) ** - alpha

def powerSpectrumHino(_freq, _z, _U, _Iu, _nonDimensional):
    """Power Spectrum (Hino)"""
    Kr = 0.0025
    m = 1.0
    n = _freq.freq

    beta = 0.0178 * alpha *  Kr * _U / _Iu ** 3.0 * (_z / 10.0) ** ((2.0 * m - 4.0) * alpha - 1.0)
    if _nonDimensional:
        Su = 0.4751 * n / beta / (1 + (n / beta) ** 2.0) ** (5.0 / 6.0)
    else:
        Su = 0.4751 / beta / (1 + (n / beta) ** 2.0) ** (5.0 / 6.0) * (_Iu * _U ) ** 2.0
    return Su

def powerSpectrumBushAndPanofsky(_freq, _z, _U, _Iw, _nonDimensional):
    """Power Spectrum (Bush & Panofsky)"""
    f_max = 0.3
    n = _freq.freq
    f = n * _z / _U
    if _nonDimensional:
        Sw = 0.6320 * f / f_max / (1.0 + 1.5 * (f / f_max) ** (5.0 / 3.0))
    else:
        Sw = 0.6320 * f / f_max / (1.0 + 1.5 * (f / f_max) ** (5.0 / 3.0)) * (_Iw * _U) ** 2.0 / n
    return Sw

def correlationFunction(_elm2, _freq):
    """Correlation Function"""
    Ruu = math.exp(- _freq.freq / _elm2.U * math.sqrt((K1 * _elm2.deltaX) ** 2.0 + (K2 * _elm2.deltaZ) ** 2.0))
    return Ruu

def aerodynamicAdmittanceDrag(_elm, _freq, _nonDimensional):
    """Aerodynamic Admittance for Drag"""
    U = _elm.U
    if _elm.type == 'girder':
        eta = _freq.freq * _elm.h / U
        XD2 = 2.0 / (K2 * eta) ** 2.0 * (K2 * eta - 1.0 + math.exp(- K2 * eta))
    elif _elm.type == 'cable':
        XD2 = 1.0

    if _nonDimensional != True:
        AD = _elm.h * _elm.l * _elm.phai
        XD2 = 4.0 * (roh * _elm.cD * AD * U ** 2.0 / 2.0) ** 2.0 * XD2 / U ** 2.0 #Ha2
    return XD2

def aerodynamicAdmittanceLift(_elm, _freq, _nonDimensional):
    """Aerodynamic Admittance for Lift"""
    a = 0.1811
    U = _elm.U
    eta = _freq.omega * _elm.b / 2.0 / U
    XL2 = (a + eta) / (a + (math.pi * a + 1.0) * eta + 2.0 * math.pi * eta ** 2.0)

    if _nonDimensional != True:
        cF = _elm.h * _elm.cD / _elm.b + _elm.dcLda
        XL2 = (roh * cF * _elm.b * _elm.l * U ** 2.0 / 2.0) ** 2.0 * XL2 / U ** 2.0
    return XL2

def mechanicalAdmittance(_omegaK, _xiSK, _xiAK, _freq):
    """Mechanical Admittance"""
    _omega = _freq.omega
    return 1.0 / ((_omegaK ** 2.0 - _omega ** 2.0) ** 2.0 + (2.0 * (_xiSK + _xiAK) * _omegaK * _omega) ** 2.0)

def aerodynamicDampingDrag(_omegaK, _modeN, _elmv):
    """Aerodynamic Damping for Drag"""
    xiAK = 0.0
    for i in _elmv:
        elm = _elmv[i]
        xiAK += elm.mode[_modeN].y ** 2.0 * roh * elm.cD * elm.h * elm.l * elm.phai * elm.U / 2.0 / _omegaK
    return xiAK

def aerodynamicDampingLift(_omegaK, _modeN, _elmv):
    """Aerodynamic Damping for Lift"""
    xiAK = 0.0
    for i in _elmv:
        elm = _elmv[i]
        cF = elm.h * elm.cD / elm.b + elm.dcLda
        xiAK += elm.mode[_modeN].z ** 2.0 * roh * cF * elm.b * elm.l * elm.U / 4.0 / _omegaK
    return xiAK

def powerSpectrumAerodynamicForce(_dof, _freq, _modeN, _elmv, _elm2m, _usePSWind):
    """Power Spectrum of Aerodynamic Force"""
    modea = []
    for i in sorted(_elmv):
    #for i in _elmv:
        elm = _elmv[i]
        if _dof == 'drag':
            modea.append(elm.mode[_modeN].y)
        elif _dof == 'lift':
            modea.append(elm.mode[_modeN].z)
    modev = numpy.matrix(modea).T
    adf = []
    ii = 0
    for i in sorted(_elmv):
    #for i in _elmv:
        elmi = _elmv[i]
        adf_temp = []
        jj = 0
        for j in sorted(_elmv):
        #for j in _elmv:
            if jj < ii:
                adf_temp.append(adf[jj][ii])
            else:
                elmj = _elmv[j]
                if _usePSWind:
                    adf_temp.append(math.sqrt(elmi.aa[_dof] * elmj.aa[_dof]) * _elm2m[i][j].corrF * _freq.psw)
                else:
                    if _dof == 'drag':
                        adf_temp.append(math.sqrt(elmi.aa[_dof] * elmj.aa[_dof]) * _elm2m[i][j].corrF * _elm2m[i][j].ps[_dof] / _freq.freq * elmi.Iu * elmi.U * elmj.Iu * elmj.U)
                    elif _dof == 'lift':
                        adf_temp.append(math.sqrt(elmi.aa[_dof] * elmj.aa[_dof]) * _elm2m[i][j].corrF * _elm2m[i][j].ps[_dof] / _freq.freq * elmi.Iu * 0.5 * elmi.U * elmj.Iu * 0.5 * elmj.U)
            jj += 1
        adf.append(adf_temp)
        ii += 1
    adfm = numpy.matrix(adf)
    psadf = modev.T * adfm * modev
    return psadf.item((0, 0))
    '''
    # Numpy is slightly faster than python list operation
    psadf = 0.0
    for i in range(len(modea)):
        for j in range(len(modea)):
            psadf += modea[i] * modea[j] * adf[i][j]

    return psadf
    '''

def powerSpectrumDisplacement(_dof, _freq, _elmv, _modev, _elm2m, _num_mode, _usePSWind):
    """Power Spectrum of Displacement"""
    #Calculate mode-dependent data

    start_time = datetime.now()
    aA_deltatime = timedelta()
    cF_deltatime = timedelta()
    pSW_deltatime = timedelta()
    mA_deltatime = timedelta()
    pSA_deltatime = timedelta()
    pSD_deltatime = timedelta()

    for i in _elmv:
        if _dof == 'drag':
            aA_start_time = datetime.now()
            _elmv[i].aa[_dof] = aerodynamicAdmittanceDrag(_elmv[i], _freq, False)
            aA_end_time = datetime.now()
            aA_deltatime += aA_end_time - aA_start_time
        elif _dof == 'lift':
            aA_start_time = datetime.now()
            _elmv[i].aa[_dof] = aerodynamicAdmittanceLift(_elmv[i], _freq, False)
            aA_end_time = datetime.now()
            aA_deltatime += aA_end_time - aA_start_time
        for j in _elmv:
            cF_start_time = datetime.now()
            _elm2m[i][j].corrF = correlationFunction(_elm2m[i][j], _freq)
            cF_end_time = datetime.now()
            cF_deltatime += cF_end_time - cF_start_time
            if _usePSWind == False:
                pSW_start_time = datetime.now()
                if _dof == 'drag':
                    _elm2m[i][j].ps[_dof] = powerSpectrumHino(_freq, _elm2m[i][j].z, _elm2m[i][j].U, _elm2m[i][j].Iu, True)
                elif _dof == 'lift':
                    _elm2m[i][j].ps[_dof] = powerSpectrumBushAndPanofsky(_freq, _elm2m[i][j].z, _elm2m[i][j].U, _elm2m[i][j].Iu * 0.5, True)
                pSW_end_time = datetime.now()
                pSW_deltatime += pS_end_time - pS_start_time

    mechAdmv = dict()
    psadfv = dict()
    psddv = dict()
    for k in range(_num_mode):
        #print(k)
        mA_start_time = datetime.now()
        omegaK = _modev[k + 1].omega
        xiSK = _modev[k + 1].damp[_dof]
        xiAK = _modev[k + 1].adamp[_dof]
        mechAdmv[k + 1] = mechanicalAdmittance(omegaK, xiSK, xiAK, _freq)
        mA_end_time = datetime.now()
        mA_deltatime += mA_end_time - mA_start_time

        pSA_start_time = datetime.now()
        psadfv[k + 1] = powerSpectrumAerodynamicForce(_dof, _freq, k + 1, _elmv, _elm2m, _usePSWind)
        pSA_end_time = datetime.now()
        pSA_deltatime += pSA_end_time - pSA_start_time

    pSD_start_time = datetime.now()
    for i in _elmv:
        elm = _elmv[i]
        psdd = 0.0
        for k in range(_num_mode):
            if _dof == 'drag':
                psdd += mechAdmv[k + 1] * psadfv[k + 1] * elm.mode[k + 1].y ** 2.0
            elif _dof == 'lift':
                psdd += mechAdmv[k + 1] * psadfv[k + 1] * elm.mode[k + 1].z ** 2.0
        psddv[i] = psdd
    pSD_end_time = datetime.now()
    pSD_deltatime += pSD_end_time - pSD_start_time

    end_time = datetime.now()
    deltatime = end_time - start_time

    print(' Total Time:' + str(deltatime.total_seconds()) + ' s')
    print(' Aerodynamic Admittance:' + str(aA_deltatime.total_seconds()) + ' s')
    print(' Correlation Function:' + str(cF_deltatime.total_seconds()) + ' s')
    print(' Power Spectrum of Wind Speed:' + str(pSW_deltatime.total_seconds()) + ' s')
    print(' Mechanical Admittance:' + str(mA_deltatime.total_seconds()) + ' s')
    print(' Power Spectrum of Aerodynamic Force:' + str(pSA_deltatime.total_seconds()) + ' s')
    print(' Power Spectrum of Displacement:' + str(pSD_deltatime.total_seconds()) + ' s')

    return psddv

def gustResponseDisplacement(_dof, _elmv, _modev, _freqv, _U10, _Iu10, _num_mode, _usePSWind, _nodePSDD):
    """Gust Response of Displacement for Drag"""
    #Calculate node-frequency-dependent data
    print('Calculation of two-point data', datetime.now().isoformat())
    elm2m = {}

    for i in _elmv:
        _elmv[i].U = Uz(_elmv[i].loc.z, _U10)
        _elmv[i].Iu = Iuz(_elmv[i].loc.z, _Iu10)
        for j in _elmv:
            if i in elm2m:
                pass
            else:
                elm2m[i] = {}
            e2 = elm2()
            e2.deltaX = _elmv[i].loc.x - _elmv[j].loc.x
            e2.deltaZ = _elmv[i].loc.z - _elmv[j].loc.z
            avZ = (_elmv[i].loc.z + _elmv[j].loc.z) / 2.0
            e2.z = avZ
            e2.U = Uz(avZ, _U10)
            e2.Iu = Iuz(avZ, _Iu10)
            elm2m[i][j] = e2

    #Calculate mode-frequency-dependent data
    print('Calculation of aerodynamic admittance', datetime.now().isoformat())
    for k in range(_num_mode):
        omegaK = _modev[k + 1].omega
        if _dof == 'drag':
            _modev[k + 1].adamp[_dof] = aerodynamicDampingDrag(omegaK, k + 1, _elmv)
        elif _dof == 'lift':
            _modev[k + 1].adamp[_dof] = aerodynamicDampingLift(omegaK, k + 1, _elmv)

    intSr = dict()
    intn2Sr = dict()
    for i in _elmv:
        intSr[i] = 0.0
        intn2Sr[i] = 0.0
    #print('calculation of power spectrum of displacement', datetime.now().isoformat(timespec = 'milliseconds'))
    print('calculation of power spectrum of displacement', datetime.now().isoformat())
    #pool = mp.Pool(2)
    for f in sorted(_freqv):
    #for f in _freqv:
        if int(f * 100) % 1 == 0:
            #print(f'{f: .4f}', datetime.now().isoformat(timespec = 'milliseconds'))
            print(f, datetime.now().isoformat())
        Sr = powerSpectrumDisplacement(_dof, _freqv[f], _elmv, _modev, elm2m, _num_mode, _usePSWind)
        for i in _elmv:
            intSr[i] += Sr[i] * _freqv[f].d_f
            intn2Sr[i] += f ** 2.0 * Sr[i] * _freqv[f].d_f
            if i == _nodePSDD:
                _freqv[f].result = Sr[i]

    print('calculation of displacement', datetime.now().isoformat())
    for i in _elmv:
        if intSr[i] != 0.0:
            sigma = math.sqrt(intSr[i])
            n = math.sqrt(intn2Sr[i]) / sigma
            root = math.sqrt(2.0 * math.log(n * T))
            G = root + 0.5772 / root
            #_elmv[i].result[_dof] = G * sigma
            _elmv[i].result[_dof] = sigma, G * sigma
        else:
            _elmv[i].result[_dof] = 0.0, 0.0
