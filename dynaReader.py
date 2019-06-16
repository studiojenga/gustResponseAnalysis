import struct, math, csv
import collections
import cProfile
import numpy
from datetime import datetime
import gustResponseAnalysis as gra
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

#print('Start analysis', datetime.now().isoformat(timespec = 'milliseconds'))
#print('Start analysis', datetime.now().isoformat())

path = 'gust_data/f13'
girder = [[8001, 8067], [8068, 8208],[8209, 8275]]
cable = [[1001, 1275], [2001, 2275]]

nodePSDD = 8138 #center of center span
#nodePSDD = 8103 # 1/4 of center span
#nodesPSDD = [8138'''center of center span''', 8103'''1/4 of center span''']

#test commit
min_f = 0.0005
#max_f = 0.2
max_f = 0.01
d_f = 0.0005

num_mode = 20
'''
min_f = 0.001
max_f = 0.05
d_f = 0.001

num_mode = 2
'''
log_damp_D = 0.03
log_damp_L = 0.03

analysisDirection = 'drag'
#analysisDirection = 'lift'

useCable = True

useWindDirection = True

usePowerSpectrumOfWind = True

class winddata:
    def __init__(self, path = '', U10 = 46.0, Iu10 = 0.1, windDir = 0.0):
        self.path = path
        self.U10 = U10
        self.Iu10 = Iu10
        self.windDir = windDir

#windData = winddata(U10 = 46.0, Iu10 = 0.13, windDir = 45.0)
#windData = winddata(U10 = 21.3876, Iu10 = 0.13)

# no angle adjusted
#windData = winddata('gust_data/W20180904125846_10min.f71', 25.9251, 0.1578 * (108.3 / 10.0) ** 0.125, 59.0487)
#windData = winddata('gust_data/W20180904130846_10min.f71', 27.1791, 0.1316 * (108.3 / 10.0) ** 0.125, 56.8334)

#windData = winddata('gust_data/W20180904131846_10min.f71', 27.509, 0.1664 * (108.3 / 10.0) ** 0.125, 54.161)
#windData = winddata('gust_data/W20180904132846_10min.f71', 26.916, 0.1516 * (108.3 / 10.0) ** 0.125, 46.363)
#windData = winddata('gust_data/W20180904135846_10min.f71', 20.540, 0.0767 * (108.3 / 10.0) ** 0.125, 11.380)

#windData = winddata('gust_data/W20180904140846_10min.f71', 23.1275, 0.0543 * (108.3 / 10.0) ** 0.125, 35.1615)
#windData = winddata('gust_data/W20180904141846_10min.f71', 22.1305, 0.0646 * (108.3 / 10.0) ** 0.125, 41.9391)
windData = winddata('gust_data/W20180904142846_10min.f71', 21.3876, 0.0562 * (108.3 / 10.0) ** 0.125, 50.7775)


class member:
    def __init__(self, _a, _b):
        self.a = _a
        self.b = _b

def readInt(_start, _data):
    return struct.unpack('I', _data[_start: _start + 4])[0]

def readFloat(_start, _data, _num = 1):
    if _num == 1:
        return struct.unpack('d', _data[_start: _start + 8])[0]
    else:
        v = []
        p = _start
        for i in range(_num):
            v.append(struct.unpack('d', _data[p: p + 8])[0])
            p += 8
        t = tuple(v)
        return t

def nodeSelector(_grd):
    selector = False
    position = 'middle'
    for i in girder:
        if _grd >= i[0] and _grd <= i[1]:
            selector = 'girder'

    if useCable:
        for i in cable:
            if _grd >= i[0] and _grd <= i[1]:
                selector = 'cable'

    return selector

#Adjust wind speed for Akashi Kaikyo Bridge
def windSpeedAdjuster(_angle):
    return math.sqrt(1.0469 / 1.9665 * math.cos(math.radians(2.368 * _angle - 33.12)) + 1.0469 / 1.9665)

if __name__ == '__main__':
    print('Start analysis', datetime.now().isoformat())

    points = dict()
    #points = collections.OrderedDict()
    members = dict()
    modes = dict()

    with open(path, 'rb') as f:
        pointer = 0
        data = f.read()
        np = readInt(pointer + 8, data)
        pointer += 132
        for i in range(np):
            grd = readInt(pointer, data)
            x, y, z = readFloat(pointer + 4, data, 3)
            xM, yM, zM, xI, yI, zI = readFloat(pointer + 32, data, 6)

            type = nodeSelector(grd)
            if type != False:
                if type == 'girder':
                    #l = 14.0
                    #h = 14.0
                    l = 14.0
                    h = 17.1
                    b = 35.5
                    phai = 0.3992
                    cD = 1.9926
                    dcLda = 1.4463
                if type == 'cable':
                    l = 14.0
                    h = 1.122
                    b = 1.122
                    phai = 1.0
                    cD = 0.700
                    dcLda = 0.0
                points[grd] = gra.elm(type = type , x = x, y = z, z = y, l = l, h = h, b = b, phai = phai, cD = cD, dcLda = dcLda)
                #points[grd] = gra.elm(type = type, loc = gra.loc(x, z, y), l = l, h = h, phai = phai, cD = cD, mode = dict(), mass = gra.dof(xM, zM, yM, xI, zI, yI))
            pointer += 164

        for i in points:
            if points[i].type == 'girder':
                for j in girder:
                    if i == j[0]:
                        l = (points[i + 1].loc.x - points[i].loc.x) * 0.5
                    elif i == j[1]:
                        l = (points[i].loc.x - points[i - 1].loc.x) * 0.5
                    elif i > j[0] and i < j[1]:
                        l = (points[i + 1].loc.x - points[i - 1].loc.x) * 0.5
            elif points[i].type == 'cable':
                for j in cable:
                    if i == j[0]:
                        l = (points[i + 1].loc.x - points[i].loc.x) * 0.5
                    elif i == j[1]:
                        l = (points[i].loc.x - points[i - 1].loc.x) * 0.5
                    elif i > j[0] and i < j[1]:
                        l = (points[i + 1].loc.x - points[i - 1].loc.x) * 0.5
            points[i].l = l

        nm = readInt(pointer + 4, data)
        pointer += 128
        for i in range(nm):
            elm = readInt(pointer, data)
            endA = readInt(pointer + 8, data)
            endB = readInt(pointer + 12, data)
            members[elm] = member(endA, endB)
            pointer += 128

        npm = readInt(pointer + 4, data)
        nmod = readInt(pointer + 8, data)
        pointer += 128

        for i in range(nmod):
            mode_number = readInt(pointer, data)
            freq, scl = readFloat(pointer + 4, data, 2)
            modes[mode_number] = gra.mode(omega = 2.0 * math.pi * freq, dampD = log_damp_D / 2.0 / math.pi, dampL = log_damp_L / 2.0 / math.pi, scale = scl)
            pointer += 164
            for j in range(npm):
                grd = readInt(pointer, data)
                x, y, z, rx, ry, rz = readFloat(pointer + 4, data, 6)
                try:
                    points[grd].mode[mode_number] = gra.dof(x, z, y, rx, rz, ry)
                except KeyError:
                    pass
                pointer += 152

        print('number of points: ', np)
        print('number of members: ', nm)
        print('number of modes: ', nmod)

    # Read power spectrum data
    if usePowerSpectrumOfWind:
        with open(windData.path, 'r') as f_ps:
            f_ps.readline()
            ps = []
            for line in f_ps:
                try:
                    ps.append([float(line[0:15]), float(line[16:30])])
                except ValueError:
                    break

        freqs = dict()
        x = []
        y = []
        y_hino = []
        i = 0
        while ps[i][0] < max_f:
            if i == 0:
                d_f = (ps[i + 1][0] - ps[i][0]) / 2.0
            elif i == len(ps) - 1:
                d_f = (ps[i][0] - ps[i - 1][0]) / 2.0
            else:
                d_f = (ps[i + 1][0] - ps[i - 1][0]) / 2.0
            if ps[i][0] != 0.0:
                freqs[ps[i][0]] = gra.freq(ps[i][0], d_f, ps[i][1])
                if useWindDirection:
                    freqs[ps[i][0]].psw = windSpeedAdjuster(windData.windDir) ** 2.0 * freqs[ps[i][0]].psw
                x.append(ps[i][0])
                y.append(ps[i][1])
                y_hino.append(gra.powerSpectrumHino(freqs[ps[i][0]], 90.0, 60.0, 0.1, False))
            i += 1

        '''
        ax = plt.gca()
        ax.loglog(x, y)
        ax.loglog(x, y_hino)
        plt.show()
        '''
    else:
        #Evenly spaced frequency
        freqs = dict()
        f = min_f
        while f < max_f:
            freqs[f] = gra.freq(f, d_f)
            f += d_f

    # Draw mode shapes
    '''
    for k in range(50):
        x = []
        y = []
        for i in points:
            elm = points[i]
            if elm.type == 'girder':
                x.append(elm.loc.x)
                y.append(elm.mode[k + 1].y)
        print('Mode: ', k + 1, ', frequency: ', modes[k + 1].omega / 2.0 / math.pi)
        #plt.cla()
        ax = plt.gca()
        ax.cla()
        ax.plot(x, y)
        ax.set_ylim(-0.03, 0.03)
        #plt.plot(x, y)
        plt.show()
    '''
    #Mechanical admittance
    '''
    for k in range(20):
        x = []
        y = []
        omegaK = modes[k + 1].omega
        dampK = modes[k + 1].damp
        add = gra.aerodynamicDampingDrag(omegaK, k + 1, points)
        for i in range(2000):
            f = i * 0.0001
            ma = gra.mechanicalAdmittance(omegaK, dampK, add, 2.0 * math.pi * f)
            x.append(f)
            y.append(ma)
        print('mode: ', k + 1, ', aerodynamic damping: ', add)
        #plt.loglog(x, y)
        plt.plot(x, y)
    plt.show()
    '''
    #Gust response analysis
    if useWindDirection:
        windData.U10 = windSpeedAdjuster(windData.windDir)* windData.U10
        #windData.Iu10 = windSpeedAdjuster(windData.windDir)* windData.Iu10

    gra.gustResponseDisplacement(analysisDirection, points, modes, freqs, windData.U10, windData.Iu10, num_mode, usePowerSpectrumOfWind, nodePSDD)

    x = []
    y = []
    for i in freqs:
        x.append(freqs[i].freq)
        y.append(freqs[i].result)
    '''
    plt.cla()
    plt.loglog(x, y)
    plt.show()
    '''

    with open('psD.csv', 'w', newline = '') as csvfile:
        csvwriter = csv.writer(csvfile)
        for i in range(len(x)):
            csvwriter.writerow([x[i], y[i]])

    x = []
    y = []
    for i in points:
        elm = points[i]
        if elm.type == 'girder':
            x.append(elm.loc.x)
            y.append(elm.result[analysisDirection])
    '''
    plt.cla()
    plt.plot(x, y)
    plt.show()
    '''
    with open('disp.csv', 'w', newline = '') as csvfile:
        csvwriter = csv.writer(csvfile)
        for i in range(len(x)):
            #csvwriter.writerow([x[i], y[i]])
            csvwriter.writerow([x[i], y[i][0], y[i][1]])

    #Equivalent mass
    '''
    for j in range(nmod):
        eqM = 0.0
        for i in points:
            eqM += points[i].mass.x * points[i].mode[j + 1].x ** 2.0
            eqM += points[i].mass.y * points[i].mode[j + 1].y ** 2.0
            eqM += points[i].mass.z * points[i].mode[j + 1].z ** 2.0
            eqM += points[i].mass.rx * points[i].mode[j + 1].rx ** 2.0
            eqM += points[i].mass.ry * points[i].mode[j + 1].ry ** 2.0
            eqM += points[i].mass.rz * points[i].mode[j + 1].rz ** 2.0

        print(j + 1, eqM)
    '''
    '''
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    for mem in members:
        x = [points[members[mem].a].loc.x, points[members[mem].b].loc.x]
        y = [points[members[mem].a].loc.y, points[members[mem].b].loc.y]
        z = [points[members[mem].a].loc.z, points[members[mem].b].loc.z]
        ax.plot(x, z, y, color = 'black')
    #plt.axes().set_aspect('equal')
    #ax.set_aspect('equal')
    #ax.plot([-2000, -2000, 4000], [-2000, 2000, 4000])
    #ax.set_xlim(1850, 2150)
    ax.set_xlim(850, 1150)
    ax.set_ylim(-150, 150)
    ax.set_zlim(0, 300)
    plt.show()
    '''
