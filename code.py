# Estudo da informação em algoritmos de leitura sensorial continua
# -*- coding: utf-8 -*-
# env: python 3.7

# TODO: 1.Media aritimetica entre os valores das informações de todos os pontos por cor (OK)
# TODO: 2.Importancia calculada entre(distancia) os pontos representativos e o novo ponto (OK)
# TODO: 3.Comparar os resultados obtidos no 1 e 2 (OK)
# TODO: 4.Efetuar um indice de disperção dos pontos de cada cor em cada eixo, e assim encontrar a relevancia de cada eixo para cada cor (OK)
# TODO: 5.Estudar cenarios previos muito bem definidos (OK)
# TODO: 6.Criação de N cenarios estimados, estimulando assim o pensamento previo ()
# TODO: 7.Usar o tempo como dimensão possivel, e levar os diferentes cenarios temporais ao surgimento de uma nova informação ()

# imports
# noinspection PyUnresolvedReferences
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
import math as mt

# functions


def euclideandist(x1, y1, z1, x2, y2, z2):
    result = mt.sqrt(mt.pow((x1 - x2), 2) + mt.pow((y1 - y2), 2) + mt.pow((z1 - z2), 2))
    return result


def sortby(elem):
    return elem[5]


def sortby2(elem):
    return elem[4]


# global variables
x, y, z, res = [], [], [], []
xg, yg, zg, xr, yr, zr, xb, yb, zb, perg, perr, perb = [], [], [], [], [], [], [], [], [], [], [], []
ax3, ax4 = 0, 0
apx, apy, apz = 0, 0, 0
matdisD, m = [], []
padd = ""

# preparing data
arr, con = [], 0
for i1 in open("data.txt", "r").readlines():
    arr += i1.replace("\n", "").split(",")
    con += 1

data = np.array(arr).reshape(con, 5)
del arr, con  # cleaning memory

for i2 in data:
    x, y, z, res = x + [i2[0]], y + [i2[1]], z + [i2[2]], res + [i2[3]]
    # calculate the representation point
    if i2[3] == "green":
        xg, yg, zg, perg = xg + [i2[0]], yg + [i2[1]], zg + [i2[2]], perg + [i2[4]]
    if i2[3] == "red":
        xr, yr, zr, perr = xr + [i2[0]], yr + [i2[1]], zr + [i2[2]], perr + [i2[4]]
    if i2[3] == "blue":
        xb, yb, zb, perb = xb + [i2[0]], yb + [i2[1]], zb + [i2[2]], perb + [i2[4]]

# printing the data frame
print("The data frame we are working with:\n", data, "\n")

# converting the string data from raw txt into float values
xg, yg, zg, perg = [float(i) for i in xg], [float(i) for i in yg], [float(i) for i in zg], [float(i) for i in perg]
xr, yr, zr, perr = [float(i) for i in xr], [float(i) for i in yr], [float(i) for i in zr], [float(i) for i in perr]
xb, yb, zb, perb = [float(i) for i in xb], [float(i) for i in yb], [float(i) for i in zb], [float(i) for i in perb]
x = [float(i) for i in x]
y = [float(i) for i in y]
z = [float(i) for i in z]

# calculating the representation position based on the weight(relevance) from each point
xmp = [np.average(xg, weights=perg)] + [np.average(xr, weights=perr)] + [np.average(xb, weights=perb)]
ymp = [np.average(yg, weights=perg)] + [np.average(yr, weights=perr)] + [np.average(yb, weights=perb)]
zmp = [np.average(zg, weights=perg)] + [np.average(zr, weights=perr)] + [np.average(zb, weights=perb)]
resmp = ['green', 'red', 'blue']

# calculation, plots configuration and op input checking
axisdisp = [0, 0.25, 0.5, 0.75, 1]
fig = ""
op = input("Add a new point? Yes(1) / No (0): ")

if op == '1':

    # config plot to new information import
    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(221, projection='3d')
    ax2 = fig.add_subplot(222, projection='3d')
    ax3 = fig.add_subplot(223, projection='3d')
    ax4 = fig.add_subplot(224, projection='3d')

    ax3.set_xlim(0, 1), ax3.set_ylim(0, 1), ax3.set_zlim(0, 1)
    ax4.set_xlim(0, 1), ax4.set_ylim(0, 1), ax4.set_zlim(0, 1)
    ax3.set_xticks(axisdisp), ax3.set_yticks(axisdisp), ax3.set_zticks(axisdisp)
    ax4.set_xticks(axisdisp), ax4.set_yticks(axisdisp), ax4.set_zticks(axisdisp)
    ax3.set_xlabel('x factor'), ax3.set_ylabel('y factor'), ax3.set_zlabel('z factor')
    ax4.set_xlabel('x factor'), ax4.set_ylabel('y factor'), ax4.set_zlabel('z factor')
    ax3.title.set_text('The New Point in data universe')
    ax4.title.set_text('The New Point in Repr. universe')

    apx = float(input("X: "))
    apy = float(input("Y: "))
    apz = float(input("Z: "))

    # 1. First Method

    print("\n1. First Method:")
    print("\nThe new information sorted by euclidean distance(each point):")
    print("Proximity -> [X, Y, Z, Universe, Point Importance, Distance, RealRelevance]")

    # calculating how important is the information for each point in universe
    matdis = []
    for h in data:
        matdis += [[h[0]] + [h[1]] + [h[2]] + [h[3]] + [h[4]] + [euclideandist(apx, apy, apz, float(h[0]), float(h[1]), float(h[2]))]]

    matdisD = sorted(matdis, key=sortby)

    for c, e in enumerate(matdisD):
        # add importance
        matdisD[c] = matdisD[c] + [((mt.sqrt(3)-e[5])/mt.sqrt(3))*float(e[4])]
        print(str(c+1) + "° -> " + str(matdisD[c]))

    relb, relr, relg = 0, 0, 0
    cb, cr, cg = 0, 0, 0

    for e in matdisD:
        if e[3] == "blue":
            relb, cb = relb + e[6], cb + 1
        if e[3] == "red":
            relr, cr = relr + e[6], cr + 1
        if e[3] == "green":
            relg, cg = relg + e[6], cg + 1

    mb, mr, mg = relb / cb, relr / cr, relg / cg
    tempm = np.array([[mb, "blue"], [mr, "red"], [mg, "green"]])
    m = tempm[np.argsort(tempm[:, 0])]
    del tempm, relb, relr, relg, cb, cr, cg, matdis
    print("\nAverage of relevance by each color:")

    for c, e in enumerate(m):
        print(str(3 - c) + "° -> " + str(m[c]))

    # 2. Second Method

    print("\n2. Second Method:")
    print("\nThe new information sorted by euclidean distance(repr point):")
    print("Proximity -> [X, Y, Z, Universe, Distance, RealRelevance]")

    # calculating how important is the information for each representation point in universe
    matdisD2, matdisD22 = [], []

    data2 = np.rot90(np.fliplr(np.array([xmp, ymp, zmp])))
    data2 = data2.tolist()

    for c, h in enumerate(data2):
        matdisD2 += [[h[0]] + [h[1]] + [h[2]] + [resmp[c]] + [euclideandist(apx, apy, apz, float(h[0]), float(h[1]), float(h[2]))]]

    matdisD22 = sorted(matdisD2, key=sortby2)

    for c, e in enumerate(matdisD22):
        # add importance
        matdisD22[c] = matdisD22[c] + [(mt.sqrt(3)-e[4])/mt.sqrt(3)]
        print(str(c+1) + "° -> " + str(matdisD22[c]))

    # 3. Axis dispersion index

    if input("\nSearch Axis dispersion by specific universe? Yes(1) / No (0): ") == '1':
        dxi, dyi, dzi = [], [], []
        for c in resmp:
            for h in data:
                if h[3] == c:
                    dxi, dyi, dzi = dxi + [["X"] + [c] + [h[0]] + [h[4]]],\
                                    dyi + [["Y"] + [c] + [h[1]] + [h[4]]],\
                                    dzi + [["Z"] + [c] + [h[2]] + [h[4]]]
        Universe = input("Which universe: (red, green, blue): ")
        print("\n3. Axis dispersion index:")
        print("\n[Axis, Universe, Position, Relevance]")
        for h in dxi:
            if Universe == h[1]:
                print(h)
        print("")
        for h in dyi:
            if Universe == h[1]:
                print(h)
        print("")
        for h in dzi:
            if Universe == h[1]:
                print(h)

    # ploting the new info
    ax3.scatter(x, y, z, zdir='z', c=res, s=15)
    ax3.scatter(apx, apy, apz, zdir='z', c='y', s=40)
    ax4.scatter(xmp, ymp, zmp, zdir='z', c=resmp, s=15)
    ax4.scatter(apx, apy, apz, zdir='z', c='y', s=40)

    for t in range(len(x)):
        ax3.plot([apx, x[t]], [apy, y[t]], [apz, z[t]], c=res[t], linestyle='dashed')
    for t in range(len(xmp)):
        ax4.plot([apx, xmp[t]], [apy, ymp[t]], [apz, zmp[t]], c=resmp[t], linestyle='dashed')

    padd = padd + str(len(data)) + "_X" + str(apx).replace(".", ",") + "-Y" + str(apy).replace(".", ",") + "-Z" + str(apz).replace(".", ",")

else:

    # config plot if no info. is added
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')

# general plot information
ax1.scatter(x, y, z, zdir='z', c=res, s=15)

ax1.title.set_text('The Data')
ax2.title.set_text('The Representation Point')
ax1.set_xlim(0, 1), ax1.set_ylim(0, 1), ax1.set_zlim(0, 1)
ax2.set_xlim(0, 1), ax2.set_ylim(0, 1), ax2.set_zlim(0, 1)
ax1.set_xticks(axisdisp), ax1.set_yticks(axisdisp), ax1.set_zticks(axisdisp)
ax2.set_xticks(axisdisp), ax2.set_yticks(axisdisp), ax2.set_zticks(axisdisp)
ax1.set_xlabel('x factor'), ax1.set_ylabel('y factor'), ax1.set_zlabel('z factor')
ax2.set_xlabel('x factor'), ax2.set_ylabel('y factor'), ax2.set_zlabel('z factor')
ax2.scatter(xmp, ymp, zmp, zdir='z', c=resmp, s=15)
'''
for r, s, t in xmp, ymp, zmp:
    ax2.text(r, s, t, "X:" + str(round(r, 2)) + "\nY:" + str(round(s, 2)) + "\nZ:" + str(round(t, 2)))'''

# free perspective plot
if input("\nOpen the free plot perspective? Yes(1) / No (0): ") == '1':
    plt.show()

# plot picture configuration
while input("\nTake a plot picture? Yes(1) / No (0): ") == '1':
    ver = input("Vertical orientation: ")
    hor = input("Horizontal orientation: ")
    ax1.view_init(elev=int(ver), azim=int(hor))
    ax2.view_init(elev=int(ver), azim=int(hor))
    if op == '1':
        ax3.view_init(elev=int(ver), azim=int(hor))
        ax4.view_init(elev=int(ver), azim=int(hor))
    plt.savefig("pictures/v" + ver + "-h" + hor + padd + ".png")

ax1.view_init(elev=30., azim=60)
ax2.view_init(elev=30., azim=60)
if op == '1':
    ax3.view_init(elev=30., azim=60)
    ax4.view_init(elev=30., azim=60)

# save data and create log file
if op == '1':
    if input("\nAdd this new info to data.txt? Yes(1) / No (0): ") == '1':
        datafile = open("data.txt", "a")
        for c in m:
            datafile.write("\n" + str(apx) + "," + str(apy) + "," + str(apz) + "," + c[1] + "," + str(c[0]))
        logfile = open("logs/log" + padd + ".txt", "w+")
        logfile.write("The new point: " + "X:" + str(apx) + " Y:" + str(apy) + " Z:" + str(apz) + "\n")
        logfile.write("The euclidean distance data:\n")
        logfile.write("Proximity -> [X, Y, Z, Color, Point Importance, Distance, RealRelevance]\n")
        for c, e in enumerate(matdisD):
            logfile.write(str(c+1) + " -> " + str(e) + "\n")
