"""Bibliothèques"""
import argparse
import math
import sys
import time
import calendar
import os

import numpy
import numpy as np

"""Constantes"""
pi=math.pi
mu= 398600.442000000
J2_earth = 0.00108262668
RT = 6378.136658
Lver = 280 #.45692
Rotearh = 360.9856474
EpVer = 36526.4992571296
RTpol = 6356.7523



def saisie(inputtxt, inputint=False, inputlen=0, inputend=0,inputdef=""):
    vi=False
    vl=False
    ve=False
    if inputint==False: vi=True
    if inputlen==0: vl = True
    if inputend==0: ve=True

    while vi==False or vl==False or ve==False:
        txt=input(inputtxt)

        if inputdef != "" and txt == "":
            txt = inputdef
            break

        if inputint == True:
            try:
                nb=int(txt)
                vi = True
            except:
                vi = False
                print(f"la saisie doit être un nombre entier")

        if inputlen != 0:
            if len(txt)<=inputlen:
                vl = True
            else:
                vl=False
                print(f"la saisie doit faire {inputlen} caractères maximum")

        if inputend != 0:
            if txt[len(txt)-len(inputend):] == inputend[:len(inputend)]:
                ve = True
            else:
                ve = False
                print(f"la saisie doit se terminer par : {inputend}")
    return txt

def readTLE():
    """Permet de lire un TLE.txt et de renvoyer les paramètres orbitaux de chaque satellite"""
    """Ouverture du fichier"""
    fichiers = os.listdir('../TLE/data/INPUT')
    for x in range(len(fichiers)): print(f"{x+1} :{fichiers[x]}",sep="\n")
    fichier = ('data/INPUT')+"/"+fichiers[int(input("Selection du fichier TLE : "))-1]

    f = open(fichier, "r")
    lignes = f.readlines()
    nb_l = int(len(lignes))

    """Variables de sortie"""
    sat = []
    Epoch = []
    INC = []
    RAAN = []
    ECC = []
    AoP = []
    MA = []
    MM = []

    """Determination du type TLE/3LE"""
    l1 =lignes[0].split()
    l2 =lignes[1].split()
    l3 =lignes[2].split()

    if nb_l%3==0:
        if l2[0]=="1" and l3[0]=="2":
            type=3
        else:
            type = 0

    elif nb_l%2==0:
        if l1[0]=="1" and l2[0]=="2":
            type=2
        else :
            type = 0
    else:
        type = 0

    """Extraction des données"""
    match type:
        case 0:
            print("merci d'utiliser un format TLE/3LE standard")
        case 2:
            print("format 2LE")
            for ligne in lignes:
                l=l+1
                if l%2==1:
                    sat[(l // 2)-1]="sat"+l//2
                    Epoch[(l // 2)-1]= "20"+ligne[6]+" "+ligne[7]
                else:
                    INC[(l // 2)-1] = ligne[2]
                    RAAN[(l // 2)-1] = ligne[3]
                    ECC[(l // 2)-1] = ligne[4]
                    AoP[(l // 2)-1] = ligne[5]
                    MA[(l // 2)-1] = ligne[6]
                    MM[(l // 2)-1] = ligne[7]
        case 3:
            print("format 3LE")
            for l in range(nb_l):
                ligne = lignes[l].split()
                if (l)%3 == 0:
                    sat.append(ligne)
                if (l) % 3 == 1:
                    #Epoch.append("20" + ligne[3][0:2] + " " +ligne[3][2:])
                    if float(ligne[3][0:2])>57:
                        yy=math.floor((100-float(ligne[3][0:2]))*365.25)
                    else: yy=math.ceil(float(ligne[3][0:2])*365.25)
                    Epoch.append(36526 + yy + float(ligne[3][2:])) #EpXLS
                    #Epoch.append(51544 + yy + float(ligne[3][2:]))  # EpXLS
                if (l) % 3 == 2:
                    INC.append(float(ligne[2]))
                    RAAN.append(float(ligne[3]))
                    ECC.append(float("0."+ligne[4]))
                    AoP.append(float(ligne[5]))
                    MA.append(float(ligne[6]))
                    MM.append(float(ligne[7]))

    #print(f"noms des satellites: {sat} \n, INC : {INC} \n, RAAN : {RAAN} \n, ECC : {ECC} \n, AoP : {AoP} \n, MA : {MA} \n, MM : {MM} \n, Epoch : {Epoch} \n, ")
    return sat, INC, RAAN, ECC, AoP, MA, MM, Epoch

def drift(ECC, INC, MM):
    period=(86400/MM)
    SMA= ((mu*period**2)/(4*pi**2))**(1/3)
    N=2*pi/period
    INC2=INC*pi/180

    precess = (-3 / 2) * (RT / SMA) * (RT / SMA) * ((N * J2_earth) / ((1 - ECC * ECC) * (1 - ECC * ECC)))

    dRAAN=(precess*math.cos(INC2))*180/pi
    dAoP =((-precess / 2) * (4 - (5 * math.sin(INC2) * math.sin(INC2))))*180/pi
    dMA = (N + (-precess / (2 * ((1 - ECC * ECC)**(1 / 2)))) * (2 - (3 * math.sin(INC2) * math.sin(INC2))))*180/pi

    return dRAAN, dAoP, dMA

def COS(x):
    a=math.cos(math.radians(x))
    return a
def SIN(x):
    a=math.sin(math.radians(x))
    return a
def TAN(x):
    a=math.tan(math.radians(x))
    return a
def ACOS(x):
    a=math.degrees(math.acos(x))
    return a
def ASIN(x):
    a=math.degrees(math.asin(x))
    return a
def ATAN(x):
    a=math.degrees(math.atan(x))
    return a
def ATAN2(x,y):
    a=math.degrees(math.atan2(x,y))
    return a

def KEP2PVT(Epoch, INC, RAAN, ECC, AoP, TA, MM):
    #convertit des coordonnées Kepleriennes en cartésiennes stellaires
    t = 86400 / MM
    a = ((t * t * mu) / (4 * pi * pi)) ** (1 / 3)
    h = ((1 - (ECC * ECC)) * a * mu) ** (1 / 2)

    p=np.array([(h * h / mu) * (1 / (1 + ECC * COS(TA))) *COS(TA),(h * h / mu) * (1 / (1 + ECC * COS(TA))) *SIN(TA),0],float)
    Vp=np.array([(mu / h) * (-SIN(TA)), (mu / h) * (ECC + COS(TA)), 0], float)

    Qx = np.array([COS(AoP) * COS(RAAN) - SIN(AoP) * COS(INC) * SIN(RAAN),-SIN(AoP) * COS(RAAN) - COS(AoP) * COS(INC) * SIN(RAAN), SIN(INC) * SIN(RAAN) ], float)
    Qy = np.array([COS(AoP) * SIN(RAAN) + SIN(AoP) * COS(INC) * COS(RAAN), -SIN(AoP) * SIN(RAAN) + COS(AoP) * COS(INC) * COS(RAAN), -SIN(INC) * COS(RAAN)],float)
    Qz =np.array([SIN(AoP) * SIN(INC), COS(AoP) * SIN(INC), COS(INC)], float)


    X = np.dot(Qx,p)
    Y = np.dot(Qy,p)
    Z = np.dot(Qz,p)
    Vx = np.dot(Qx,Vp)
    Vy = np.dot(Qy,Vp)
    Vz = np.dot(Qz,Vp)

    return Epoch, X, Y, Z, Vx, Vy, Vz

def PVT2LLA(Epoch, X, Y, Z, Vx, Vy, Vz):
    Rayon = math.sqrt(X**2+Y**2+Z**2)
    Lat = ASIN(Z/Rayon)
    LongT = (ATAN2(Y, X)-(Lver+Rotearh*(Epoch-EpVer)))%360         #=Mod_ang_180°(Ra - (LonVernal1Jan + (AngEarthRot24h * (epoch - J2000))))    (angle%360)-360*((angle%360)//180)
    Long = LongT-360*(LongT//180)
    Alt = Rayon-RT+(RT-RTpol)*SIN(Lat)**2
    return Epoch, Lat, Long, Alt #, bearing, slope

def propa_TLE(Ep0, occ, incr, INC, RAAN, ECC, AoP, MA, MM, EpochTLE):
    KEPSAT=[]
    for x in range(len(INC)):
        dRAAN, dAoP, dMA = drift(ECC[x], INC[x], MM[x])
        KEPT = []
        INC2=INC[x]
        ECC2=ECC[x]
        MM2=MM[x]
        z=(Ep0-EpochTLE[x])*86400/incr
        for y in range(occ+1):
            RAAN2=RAAN[x]+dRAAN * incr * (y+z)
            AoP2 = AoP[x] + dAoP * incr * (y+z)
            MA2 = (MA[x] + dMA * incr * (y+z))
            Ep2 = Ep0 + (incr * (y))/86400
            KEPT.append((Ep2, INC2, RAAN2, ECC2, AoP2, MA2, MM2))
        KEPSAT.append(KEPT)
        #print(*KEPSAT[x], sep="\n")
    return KEPSAT


def menu():
    NME, INC, RAAN, ECC, AoP, MA, MM, EpochTLE = readTLE()
    Ep0 = 25569 + calendar.timegm(time.strptime(input("Date de début de l'étude au format JJ/MM/YYYY HH:MM:SS : "),'%d/%m/%Y %H:%M:%S'))/86400
    incr= int(saisie(f"Incrément en secondes : ",inputint=True))
    Ep1 = 25569 + calendar.timegm(time.strptime(input("Date de fin de l'étude au format JJ/MM/YYYY HH:MM:SS : "), '%d/%m/%Y %H:%M:%S')) / 86400
    occ= int(((Ep1-Ep0)*86400)/incr)

    KEPSAT = propa_TLE(Ep0,  occ, incr, INC, RAAN, ECC, AoP, MA, MM, EpochTLE)
    LLA = False
    if input("Format de sortie ? PVT ou LLA : " ) == "LLA": LLA = True
    #satoutput= []
    for x in range(len(KEPSAT)):
        print(NME[x],INC[x],RAAN[x], ECC[x], AoP[x], MA[x], MM[x], EpochTLE[x])
        occoutput = []
        for y in range(len(KEPSAT[0])):
                PVT=KEP2PVT(*KEPSAT[x][y])
                if LLA == True:
                    occoutput.append(PVT2LLA(*PVT))
                else:
                    occoutput.append(PVT)
        #satoutput.append(occoutput)
        #print(occoutput[0])

        fichier = "../TLE/data/OUTPUT/OEM_"+str(NME[x])+".txt"
        f = open(fichier, "w")
        f.write(f"CIC_OEM_VERS \t= 2.0\nCREATION_DATE  \t= {time.strftime('%Y-%m-%dT%H:%M:%SZ',time.localtime())}\nORIGINATOR     \t= e.NOVA\n\nMETA_START\n\nOBJECT_NAME \t= {NME[x]}\nOBJECT_ID \t= {NME[x]}\n\nCENTER_NAME \t= EARTH\nREF_FRAME   \t= EME2000\nTIME_SYSTEM \t= UTC\n\nMETA_STOP\n\n")
        for y in range(len(KEPSAT[0])):
            if LLA == False : txt=f"{int(occoutput[y][0])+15018} {round((occoutput[y][0]-int(occoutput[y][0]))*86400,2)}\t{occoutput[y][1]} {occoutput[y][2]} {occoutput[y][3]} {occoutput[y][4]} {occoutput[y][5]} {occoutput[y][6]}\n"
            else: txt = f"{int(occoutput[y][0])+15018} {round((occoutput[y][0] - int(occoutput[y][0])) * 86400, 2)}\t{occoutput[y][1]} {occoutput[y][2]} {occoutput[y][3]}\n"
            f.write(txt)
        f.close()
    print(time.strftime("%d/%m/%Y %Hh%Mm%Ss",time.localtime()))



#main
menu()
