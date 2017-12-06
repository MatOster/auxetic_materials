def laenge(vec):
    l = 0
    for i in range(len(vec)):
        l+=vec[i]**2
    #l=np.sqrt(l)
    return l


#determing the energy and its first derivative


from sympy import Symbol
import numpy as np
import sympy
from sympy.tensor.array import derive_by_array

delta = 1

x3 = Symbol("x3")
x4 = Symbol("x4")
x5 = Symbol("x5")
x6 = Symbol("x6")

x8 = Symbol("x8")

x12 = Symbol("x12")
x13 = Symbol("x13")
x14 = Symbol("x14")

x18 = Symbol("x18")
x19 = Symbol("x19")
x20 = Symbol("x20")

x24 = Symbol("x24")
x25 = Symbol("x25")
x26 = Symbol("x26")
x27 = Symbol("x27")
x28 = Symbol("x28")
x29 = Symbol("x29")

x36 = Symbol("x36")
x37 = Symbol("x37")
x38 = Symbol("x38")
x39 = Symbol("x39")
x40 = Symbol("x40")
x41 = Symbol("x41")
x42 = Symbol("x42")
x43 = Symbol("x43")
x44 = Symbol("x44")
x45 = Symbol("x45")
x46 = Symbol("x46")
x47 = Symbol("x47")

x57 = Symbol("x57")
x58 = Symbol("x58")
x59 = Symbol("x59")

x66 = Symbol("x66")
x67 = Symbol("x67")
x68 = Symbol("x68")

x72 = Symbol("x72")
x73 = Symbol("x73")
x74 = Symbol("x74")

x81 = Symbol("x81")
x82 = Symbol("x82")
x83 = Symbol("x83")

x93 = Symbol("x93")
x94 = Symbol("x94")
x95 = Symbol("x95")

x102 = Symbol("x102")
x103 = Symbol("x103")
x104 = Symbol("x104")

x108 = Symbol("x108")
x109 = Symbol("x109")
x110 = Symbol("x110")

x114 = Symbol("x114")
x115 = Symbol("x115")
x116 = Symbol("x116")

x123 = Symbol("x123")
x124 = Symbol("x124")
x125 = Symbol("x125")

x138 = Symbol("x138")
x139 = Symbol("x139")
x140 = Symbol("x140")

x144 = Symbol("x144")
x145 = Symbol("x145")
x146 = Symbol("x146")

x159 = Symbol("x159")
x160 = Symbol("x160")
x161 = Symbol("x161")

#Periodicity in z and y direction
gam = Symbol("gam")
epsilon = Symbol("epsilon")

#Lagrange Multiplicators
lambda1=Symbol("lambda1")
lambda2=Symbol("lambda2")
lambda3=Symbol("lambda3")
lambda4=Symbol("lambda4")
lambda5=Symbol("lambda5")
lambda6=Symbol("lambda6")
lambda7=Symbol("lambda7")
lambda8=Symbol("lambda8")
lambda9=Symbol("lambda9")
lambda10=Symbol("lambda10")
lambda11=Symbol("lambda11")
lambda12=Symbol("lambda12")
lambda13=Symbol("lambda13")
lambda14=Symbol("lambda14")
lambda15=Symbol("lambda15")
lambda16=Symbol("lambda16")
lambda17=Symbol("lambda17")
lambda18=Symbol("lambda18")
lambda19=Symbol("lambda19")
lambda20=Symbol("lambda20")
lambda21=Symbol("lambda21")
lambda22=Symbol("lambda22")
lambda23=Symbol("lambda23")
lambda24=Symbol("lambda24")
lambda25=Symbol("lambda25")
lambda26=Symbol("lambda26")
lambda27=Symbol("lambda27")
lambda28=Symbol("lambda28")
lambda29=Symbol("lambda29")
lambda30=Symbol("lambda30")
lambda31=Symbol("lambda31")
lambda32=Symbol("lambda32")
lambda33=Symbol("lambda33")
lambda34=Symbol("lambda34")
lambda35=Symbol("lambda35")
lambda36=Symbol("lambda36")
lambda37=Symbol("lambda37")
lambda38=Symbol("lambda38")
lambda39=Symbol("lambda39")
lambda40=Symbol("lambda40")
lambda41=Symbol("lambda41")
lambda42=Symbol("lambda42")


x0 = 0.9
x1 = 0.5
x2 = 0.75
x21 = 0.35
x22 = 0.5
x23 = 0.75
x7 = 0.5

#Energy to be minimized
g0 = (x0 - x3)**2 + (x1 - x4)**2 + (-x5-gam + x2)**2+ (x0 - x12)**2 + (x1 - x13)**2 \
     + (-x14 + x2)**2+ (x21 - x45)**2 + (x22 - x46)**2 + (x23 - x47)**2+ (x21 - x36)**2 \
     + (x22 - x37)**2 + (x23 - x38)**2+(x72 - x93)**2 + (x73 - x94)**2 + (x74 - x95)**2+(-x123 + x45)**2 \
     + (-x124 + x46)**2 + (-x125 + x47)**2+(-x123 + x36)**2 + (-x124-epsilon + x37)**2 \
     + (-x125 + x38)**2+(x93 - x81)**2 + (x94+epsilon - x82)**2 + (x95 - x83-gam)**2+(x18 - x42)**2 + \
     (x19 - x43)**2 + (x20 - x44)**2+(-x114 + x42)**2 + (-x115 + x43)**2 + (-x116 + x44)**2\
     +(x18 - x159-delta)**2 + (x19 - x160)**2 + (x20 - x161)**2+(-x12 + x6)**2 + (-x13 + x7)**2 + (-x14 + x8)**2\
     +(x3 - x6)**2 + (x4 - x7)**2 + (x5 - x8)**2+(x24 - x57)**2 + (x25 - x58)**2 + (x26 - x59)**2+(x27 - x39)**2 \
     + (x28 - x40)**2 + (x29 - x41)**2+(x114 - x159)**2 + (x115 - x160)**2 + (x116 - x161)**2+(x108 - x39+delta)**2 \
     + (x109+epsilon - x40)**2 + (x110 - x41)**2+(x66 - x72)**2 + (x67 - x73)**2 + (x68 - x74)**2+(-x108 + x138+delta)**2 \
     + (-x109 + x139)**2 + (-x110 + x140)**2+(-x144+delta + x24)**2 + (-x145 + x25)**2 + (-x146 + x26+gam)**2+(x102 - x144)**2 \
     + (x103 - x145)**2 + (x104+gam - x146)**2+(x138 - x27)**2 + (x139 - x28)**2 + (x140 - x29)**2+(x81 - x66)**2 \
     + (x82 - x67-epsilon)**2 + (x83 - x68)**2+(-x102 + x57)**2 + (-x103 + x58)**2 + (-x104 + x59)**2 +(x0-x21)**2+(x1-x22)**2\
     +(x2-x23)**2 +(x93-x123)**2+(x94-x124)**2+(x95-x125)**2+ (x36-x42)**2+(x37-x43)**2+(x38-x44)**2\
     +(x12-x18)**2+(x13-x19)**2+(x14-x20)**2+(x3-x24)**2+(x4-x25)**2+(x5-x36)**2+(x39-x159)**2+(x40-x160)**2\
     +(x41-x161)**2+(x72-x114)**2+(x73+epsilon-x115)**2+(x74-x116)**2+(x66+delta-x108)**2+(x67-x109)**2+(x68-x110)**2\
     +(x144-x138)**2+(x145-x139)**2+(x146-x140-gam)**2+(x81-x102)**2+(x82-x103-epsilon)**2+(x83-x104)**2+(x57-x45)**2\
     +(x58-x46)**2+(x59+gam-x47)**2+(x6-x27-delta)**2+(x7-x28)**2+(x8-x29)**2

#constraints
K=[((x93-x123)**2+(x94-x124)**2+(x95-x125)**2-0.3025),
   ((x36-x42)**2+(x37-x43)**2+(x38-x44)**2-0.3025),
   ((x12-x18)**2+(x13-x19)**2+(x14-x20)**2-0.3025),
   ((x3-x24)**2+(x4-x25)**2+(x5-x36)**2-0.3025),
   ((x39-x159)**2+(x40-x160)**2+(x41-x161)**2-0.3025),
   ((x72-x114)**2+(x73+epsilon-x115)**2+(x74-x116)**2-0.3025),
   ((x66+delta-x108)**2+(x67-x109)**2+(x68-x110)**2-0.3025),
   ((x144-x138)**2+(x145-x139)**2+(x146-x140-gam)**2-0.3025),
   ((x81-x102)**2+(x82-x103-epsilon)**2+(x83-x104)**2-0.3025),
   ((x57-x45)**2+(x58-x46)**2+(x59+gam-x47)**2-0.3025),
   ((x6-x27-delta)**2+(x7-x28)**2+(x8-x29)**2-0.3025),
   (x1*(x14 + x5+gam) - x2*(x13 + x4) + x22*(x38 + x47) - x23*(x37 + x46)),
   (-x0*(x14 + x5+gam) + x2*(x12 + x3) - x21*(x38 + x47) + x23*(x36 + x45)),
   (x0*(x13 + x4) - x1*(x12 + x3) + x21*(x37 + x46) - x22*(x36 + x45)),
   (x124*(x38 + x47) - x125*(x37-epsilon + x46) + x94*(x74 + x83+gam) - x95*(x73 + x82-epsilon)),
   (-x123*(x38 + x47) + x125*(x36 + x45) - x93*(x74 + x83+gam) + x95*(x72 + x81)),
   (x123*(x37-epsilon + x46) - x124*(x36 + x45) + x93*(x73 + x82-epsilon) - x94*(x72 + x81)),
   (x37*(x125 + x23) - x38*(x124+epsilon + x22) + x43*(x116 + x20) - x44*(x115 + x19)),
   (-x36*(x125 + x23) + x38*(x123 + x21) - x42*(x116 + x20) + x44*(x114 + x18)),
   (x36*(x124+epsilon + x22) - x37*(x123 + x21) + x42*(x115 + x19) - x43*(x114+ x18)),
   (x13*(x2 + x8) - x14*(x1 + x7) + x19*(x161 + x44) - x20*(x160 + x43)),
   (-x12*(x2 + x8) + x14*(x0 + x6) - x18*(x161 + x44) + x20*(x159+delta + x42)),
   (x12*(x1 + x7) - x13*(x0 + x6) + x18*(x160 + x43) - x19*(x159+delta + x42)),
   (x25*(x146 + x59) - x26*(x145 + x58) + x4*(x2-gam + x8) - x5*(x1 + x7)),
   (-x24*(x146 + x59) + x26*(x144 + x57) - x3*(x2-gam + x8) + x5*(x0 + x6)),
   (x24*(x145 + x58) - x25*(x144 + x57) + x3*(x1 + x7) - x4*(x0 + x6)),
   (x160*(x116 + x20) - x161*(x115 + x19) + x40*(x110 + x29) - x41*(x109+epsilon + x28)),
   (-x159*(x116 + x20) + x161*(x114 + x18) - x39*(x110 + x29) + x41*(x108-delta + x27)),
   (x159*(x115 + x19) - x160*(x114 + x18) + x39*(x109+epsilon + x28) - x40*(x108-delta + x27)),
   (x115*(x161 + x44) - x116*(x160 + x43) + x73*(x68+gam + x95+gam) - (x74+gam)*(x67 + x94)),
   (-x114*(x161 + x44) + x116*(x159 + x42) - x72*(x68+gam + x95+gam) + (x74+gam)*(x66 + x93)),
   (x114*(x160 + x43) - x115*(x159 + x42) + x72*(x67 + x94) - x73*(x66 + x93)),
   (x109*(x41 + x140) - x110*(x40 + x139) + x67*(x74 + x83) - x68*(x73 + x82-epsilon)),
   (-x108*(x41 + x140) + x110*(x39 + x138+delta) - (x66+delta)*(x74 + x83) + x68*(x72+delta + x81+delta)),
   (x108*(x40 + x139) - x109*(x39 + x138+delta) + (x66+delta)*(x73 + x82-epsilon) - x67*(x72+delta + x81+delta)),
   (x139*(x110+gam + x29+gam) - (x140+gam)*(x109 + x28) + x145*(x104+gam + x26+gam) - x146*(x103 + x25)),
   (-x138*(x110+gam + x29+gam) + x140*(x108-delta + x27) - x144*(x104+gam + x26+gam) + x146*(x102 + x24-delta)),
   (x138*(x109 + x28) - x139*(x108-delta + x27) + x144*(x103 + x25) - x145*(x102 + x24-delta)),
   (x103*(x146 + x59) - x104*(x145 + x58) + x82*(x68 + x95-gam) - x83*(x67+epsilon + x94+epsilon)),
   (-x102*(x146 + x59) + x104*(x144 + x57) - x81*(x68 + x95-gam) + x83*(x66 + x93)),
   (x102*(x145 + x58) - x103*(x144 + x57) + x81*(x67+epsilon + x94+epsilon) - x82*(x66 + x93)),
   (x58*(x104 + x26) - (x59+gam)*(x103 + x25)),
   (-x57*(x104 + x26) + x59*(x102 + x24)),
   (x57*(x103 + x25) - x58*(x102 + x24)),
   (x28*(x140 + x41) - x29*(x139 + x40) + x7*(x14 + x5) - x8*(x13 + x4)),
   (-(x27+delta)*(x140 + x41) + x29*(x138+delta + x39+delta) - x6*(x14 + x5) + x8*(x12 + x3)),
   ((x27+delta)*(x139 + x40) - x28*(x138+delta + x39+delta) + x6*(x13 + x4) - x7*(x12 + x3))]


lambdaall = [lambda1,
lambda2,
lambda3,
lambda4,
lambda5,
lambda6,
lambda7,
lambda8,
lambda9,
lambda10,
lambda11,
lambda12,
lambda13,
lambda14,
lambda15,
lambda16,
lambda17,
lambda18,
lambda19,
lambda20,
lambda21,
lambda22,
lambda23,
lambda24,
lambda25,
lambda26,
lambda27,
lambda28,
lambda29,
lambda30,
lambda31,
lambda32,
lambda33,
lambda34,
lambda35,
lambda36,
lambda37,
lambda38,
lambda39,
lambda40,
lambda41,
lambda42]

g1 = 0

for i in range(42):
 g1 += K[i]*lambdaall[i]

lagrange = g0+g1

variables = [
x3,
x4, 
x5, 
x6,
x8,
x12,
x13,
x14,
x18,
x19,
x20,
x24,
x25,
x26,
x27,
x28,
x29,
x36,
x37,
x38,
x39,
x40,
x41,
x42,
x43,
x44,
x45,
x46,
x47,
x57,
x58,
x59,
x66,
x67,
x68,
x72,
x73,
x74,
x81,
x82,
x83,
x93,
x94,
x95,
x102,
x103,
x104,
x108,
x109,
x110,
x114,
x115,
x116,
x123,
x124,
x125,
x138,
x139,
x140,
x144,
x145,
x146,
x159,
x160,
x161,
epsilon,
gam,
lambda1,
lambda2,
lambda3,
lambda4,
lambda5,
lambda6,
lambda7,
lambda8,
lambda9,
lambda10,
lambda11,
lambda12,
lambda13,
lambda14,
lambda15,
lambda16,
lambda17,
lambda18,
lambda19,
lambda20,
lambda21,
lambda22,
lambda23,
lambda24,
lambda25,
lambda26,
lambda27,
lambda28,
lambda29,
lambda30,
lambda31,
lambda32,
lambda33,
lambda34,
lambda35,
lambda36,
lambda37,
lambda38,
lambda39,
lambda40,
lambda41,
lambda42
]


# optimizing under constraints by using gradient descent on lagrange function

# initial value
xstart= [
    0.75, 0.65, 0.0,
    0.6, 0.25,
    0.75, 0.35, 0.5,
    0.75, 0.9, 0.5,
    0.75, 0.1, 0.0,
    0.15, 0.5, 0.25,
    0.5, 0.75, 0.9,
    0.0, 0.75, 0.1,
    0.5, 0.75, 0.35,
    0.5, 0.25, 0.6,
    0.5, 0.25, 0.15,
    0.4, 0.0, 0.25,
    0.25, 0.15, 0.5,
    0.25, 0.85, 0.0,
    0.1, 0.0, 0.75,
    0.25, 0.4, 0.0,
    0.85, 0.0, 0.25,
    0.25, 0.6, 0.5,
    0.65, 0.0, 0.75,
    0.0, 0.25, 0.4,
    0.0, 0.25, 0.85,
    0.0, 0.75, 0.65,
    1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1, 1, 1, 1,
    1, 1, 1
]

#gradient descent
def grad_des(f,x, var):
    i = 0
    eta=1
    max_It = 10000
    x=np.array(x)
    deriv = derive_by_array(f, var )
    while(i<max_It):# and np.linalg.norm(f(x))>10**(-10)):
        grad = np.array(deriv.subs([(x3,x[0]),
                        (x4,x[1]),
                        (x5,x[2]),
                        (x6,x[3]),
                        (x8,x[4]),
                        (x12,x[5]),
                        (x13,x[6]),
                        (x14,x[7]),
                        (x18,x[8]),
                        (x19,x[9]),
                        (x20,x[10]),
                        (x24,x[11]),
                        (x25,x[12]),
                        (x26,x[13]),
                        (x27,x[14]),
                        (x28,x[15]),
                        (x29,x[16]),
                        (x36,x[17]),
                        (x37,x[18]),
                        (x38,x[19]),
                        (x39,x[20]),
                        (x40,x[21]),
                        (x41,x[22]),
                        (x42,x[23]),
                        (x43,x[24]),
                        (x44,x[25]),
                        (x45,x[26]),
                        (x46,x[27]),
                        (x47,x[28]),
                        (x57,x[29]),
                        (x58,x[30]),
                        (x59,x[31]),
                        (x66,x[32]),
                        (x67,x[33]),
                        (x68,x[34]),
                        (x72,x[35]),
                        (x73,x[36]),
                        (x74,x[37]),
                        (x81,x[38]),
                        (x82,x[39]),
                        (x83,x[40]),
                        (x93,x[41]),
                        (x94,x[42]),
                        (x95,x[43]),
                        (x102,x[44]),
                        (x103,x[45]),
                        (x104,x[46]),
                        (x108,x[47]),
                        (x109,x[48]),
                        (x110,x[49]),
                        (x114,x[50]),
                        (x115,x[51]),
                        (x116,x[52]),
                        (x123,x[53]),
                        (x124,x[54]),
                        (x125,x[55]),
                        (x138,x[56]),
                        (x139,x[57]),
                        (x140,x[58]),
                        (x144,x[59]),
                        (x145,x[60]),
                        (x146,x[61]),
                        (x159,x[62]),
                        (x160,x[63]),
                        (x161,x[64]),
                        (epsilon,x[65]),
                        (gam,x[66]),
                        (lambda1,x[67]),
                        (lambda2,x[68]),
                        (lambda3,x[69]),
                        (lambda4,x[70]),
                        (lambda5,x[71]),
                        (lambda6,x[72]),
                        (lambda7,x[73]),
                        (lambda8,x[74]),
                        (lambda9,x[75]),
                        (lambda10,x[76]),
                        (lambda11,x[77]),
                        (lambda12,x[78]),
                        (lambda13,x[79]),
                        (lambda14,x[80]),
                        (lambda15,x[81]),
                        (lambda16,x[82]),
                        (lambda17,x[83]),
                        (lambda18,x[84]),
                        (lambda19,x[85]),
                        (lambda20,x[86]),
                        (lambda21,x[87]),
                        (lambda22,x[88]),
                        (lambda23,x[89]),
                        (lambda24,x[90]),
                        (lambda25,x[91]),
                        (lambda26,x[92]),
                        (lambda27,x[93]),
                        (lambda28,x[94]),
                        (lambda29,x[95]),
                        (lambda30,x[96]),
                        (lambda31,x[97]),
                        (lambda32,x[98]),
                        (lambda33,x[99]),
                        (lambda34,x[100]),
                        (lambda35,x[101]),
                        (lambda36,x[102]),
                        (lambda37,x[103]),
                        (lambda38,x[104]),
                        (lambda39,x[105]),
                        (lambda40,x[106]),
                        (lambda41,x[107]),
                        (lambda42,x[108])]))
        print(list(grad))
        print(laenge(list(grad)))
        x=x-0.005*grad

        #print(len(grad))
        #print(len(x))
        i+=1
    return x

print(grad_des(lagrange,np.array(xstart),variables))