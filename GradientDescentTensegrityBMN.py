def laenge(vec):
    l = 0
    for i in range(len(vec)):
        l += vec[i] ** 2
    # l=np.sqrt(l)
    return l


# determing the energy and its first derivative
import numpy as np
import Function_Tensegrity_BMN

delta = 1
precision = 10**(-20)


x0 = 0.9
x1 = 0.5
x2 = 0.75
x21 = 0.35
x22 = 0.5
x23 = 0.75
x7 = 0.5








# optimizing under constraints by using gradient descent on lagrange function

# initial value
xstart = [0.6350662, 0.34146648, 0.52486533, 0.81786748, 0.66528335, 0.82054711,
          0.47329207, 0.76864461, 1.08040507, 0.6138582, 0.30430362, 0.15579433,
          0.11530435, 0.09596973, 0.20286325, 0.28529558, 0.34124068, 0.37463028,
          0.5329729, 0.72472759, 0.74810181, 0.34799515, 0.16167056, 0.77639817,
          0.5940497, 0.35564275, 0.25281303, 0.35498424, 0.61267202, 0.23869524,
          0.20072191, 0.16191315, -0.04380109, 0.23039896, 0.27362851, -0.0223216,
          0.28391583, 0.29417746, -0.07459657, 0.4036394, 0.3323921, -0.03591548,
          0.27561647, 0.27536401, 0.32646284, 0.03227997, 0.04326033, 0.43133938,
          0.16568436, 0.12045672, 0.50459507, 0.56113839, 0.35321082, 0.31530871,
          0.38933585, 0.68419579, -0.10309803, 0.05607774, 0.11948611, 0.42729827,
          -0.07383446, -0.01311121, 0.29527806, 0.57403447, 0.37746696, 0.13020389,
          -0.07688131, -1.05019491, -0.99254586, -1.07530414, 0.37814361, -0.87703227,
          -0.8103113, -0.85205048, 0.04922497, -1.02207783, -0.37453097, -0.02102562,
          0.27347649, 0.94397776, 0.80228034, -0.17010101, -0.05496047, -0.0837463,
          -0.01010334, 0.29465981, 0.36552965, 0.09929671, 0.3712475, 0.33290934,
          0.4918331, 0.97489573, 0.72193499, -0.5010487, -0.16138041, 0.14591098,
          -0.34587609, 0.09419332, 0.21411847, 0.94060483, 0.40992572, 0.27225244,
          0.2381762, -1.24691014, 0.89592576, -0.32270434, 1.11581378, 0.99298145,
          -5.1490275]


# gradient descent
def grad_des(x):
    i = 0
    max_It = 1009
    x = np.array(x)
    xnow = 0 * x
    gradnow = 0 * x
    norm =1
    while (i < max_It and norm > precision:  
        grad = Function_Tensegrity_BMN.gradient(x,delta)
        norm =laenge(list(grad))
        xprev = xnow
        xnow = x
        s = xnow - xprev
        gradprev = gradnow
        gradnow = grad
        g = gradnow - gradprev
        alpha1 = .01 * s.dot(g) / laenge(g)
        alpha2 = .01 * laenge(s) / s.dot(g)
        if i > 10: alpha1 = .5 * (alpha1 + alpha2)
        x = x - alpha1 * grad
        i += 1
    return x

print(grad_des(np.array(xstart))
