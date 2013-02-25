from __future__ import with_statement
import pdb;
# Adapted from matlab implementation...
#
#OKADA85 Surface deformation due to a finite rectangular source.
#    [uN,uE,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(E,N,/
#        DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN,NU)
#    computes Okada's 1985 solution for displacements, tilts and strains in a
#    geographic referential (East,North,Up). The fault top edge is centered
#    at coordinates (0,0).
#
#        E, N    : Matrix coordinates of observation points
#        DEPTH   : Depth of the fault top edge (DEP > 0)
#        STRIKE  : Strike-angle from North (in degrees)
#        DIP     : Dip-angle (in degrees)
#        LENGTH  : Fault length in strike direction (LEN > 0)
#        WIDTH   : Fault width in dip direction (WIDTH > 0)
#        RAKE    : Slip-angle direction on fault plane (in degrees)
#        SLIP    : Dislocation in rake direction
#        OPEN    : Dislocation in tensile component
#        NU      : Poisson's ratio
#
#    returns the following variables (same size as E and N):
#    uN,uE,uZ        : Displacements (Unit of SLIP and OPEN)
#    uZE,uZN         : Tilts (in radian)
#    uNN,uNE,uEN,uEE : Strains (Unit of SLIP and OPEN)/(Unit of N,E,..,WIDTH)
#
#    It is also possible to produce partial outputs, with following syntax:
#        [uN,uE,uZ] = OKADA85(...) for displacements only
#        [uN,uE,uZ,uZE,uZN] = OKADA85(...) for displacements and tilts
#        [uN,uE,uZ,uNN,uNE,uEN,uEE] = OKADA85(...) for displacements and strains
#        [uZE,uZN] = OKADA85(...) for tilts only
#        [uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(...) for tilts and strains
#        [uNN,uNE,uEN,uEE] = OKADA85(...) for strains only
#
#    Formulas and notations similar to Okada [1985]. Convention for fault
#    parameters from Aki & Richards [1980].
#
#    Author: Francois Beauducel <beauducelipgp.fr>
#      Institut de Physique du Globe de Paris, 1997-2009.
#
#    References:
#    Okada Y., Surface deformation due to shear and tensile faults in a
#       half-space, Bull. Seismol. Soc. Am., 75:4, 1135-1154, 1985.
#    Aki K., and P. G. Richards, Quantitative seismology, Freemann & Co,
#      New York, 1980.
    
    
#    Copyright (c) 1997-2009, Francois Beauducel, covered by BSD License.
#    All rights reserved.
#
#    Redistribution and use in source and binary forms, with or without 
#    modification, are permitted provided that the following conditions are 
#    met:
#
#       * Redistributions of source code must retain the above copyright 
#         notice, this list of conditions and the following disclaimer.
#       * Redistributions in binary form must reproduce the above copyright 
#         notice, this list of conditions and the following disclaimer in 
#         the documentation and/or other materials provided with the distribution
#                               
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#    POSSIBILITY OF SUCH DAMAGE.
    
from math import sqrt, sin, cos, atan, log, pi

eps=0.0000001
    
def okada85(e,n,depth,strike,dip,L,W,U1,U2,U3,nu=0.25):
    '''
    Calculate the okada displacement.  Input is

    e      East coordinate
    n      North coordinate
    depth  Depth of fault reference point
    strike Strike of fault trace, measuring north from east
    dip    Dip of fault, measuring up to left along fault
    length Length of fault along strike
    width  Width of fault down dip
    U1     Strike slip component
    U2     Dip slip component
    U3     Tensile (opening) component
    nu     Poissons ratio, default 0.25
    '''
    # From degrees to radians
    strike = strike*pi/180
    delta = dip*pi/180
    
    # Converts fault coordinates into Okada's reference system
    x = sin(strike)*n + cos(strike)*e
    y = cos(strike)*n - sin(strike)*e
    
    # Variable substitution (independent from xi and eta)
    p = y*cos(delta) + depth*sin(delta)
    q = y*sin(delta) - depth*cos(delta)
    
    # Displacements
    ux = -U1/(2*pi) * chinnery(ux_ss,x,p,L,W,q,delta,nu) \
        - U2/(2*pi) * chinnery(ux_ds,x,p,L,W,q,delta,nu) \
        + U3/(2*pi) * chinnery(ux_tf,x,p,L,W,q,delta,nu)
    
    uy = -U1/(2*pi) * chinnery(uy_ss,x,p,L,W,q,delta,nu) \
        - U2/(2*pi) * chinnery(uy_ds,x,p,L,W,q,delta,nu) \
        + U3/(2*pi) * chinnery(uy_tf,x,p,L,W,q,delta,nu)
    
    uz = -U1/(2*pi) * chinnery(uz_ss,x,p,L,W,q,delta,nu) \
        - U2/(2*pi) * chinnery(uz_ds,x,p,L,W,q,delta,nu) \
        + U3/(2*pi) * chinnery(uz_tf,x,p,L,W,q,delta,nu)
    
    # Rotation from Okada's axes to geographic
    ue = cos(strike)*ux - sin(strike)*uy
    un = sin(strike)*ux + cos(strike)*uy
    
    # Tilt
    uzx = -U1/(2*pi) * chinnery(uzx_ss,x,p,L,W,q,delta,nu) \
         - U2/(2*pi) * chinnery(uzx_ds,x,p,L,W,q,delta,nu) \
         + U3/(2*pi) * chinnery(uzx_tf,x,p,L,W,q,delta,nu)
    
    uzy = -U1/(2*pi) * chinnery(uzy_ss,x,p,L,W,q,delta,nu) \
         - U2/(2*pi) * chinnery(uzy_ds,x,p,L,W,q,delta,nu) \
         + U3/(2*pi) * chinnery(uzy_tf,x,p,L,W,q,delta,nu)
    
    # Rotation from Okada's axes to geographic
    uze = cos(strike)*uzx - sin(strike)*uzy
    uzn = sin(strike)*uzx + cos(strike)*uzy
    
    # Strain
    uxx = -U1/(2*pi) * chinnery(uxx_ss,x,p,L,W,q,delta,nu) \
         - U2/(2*pi) * chinnery(uxx_ds,x,p,L,W,q,delta,nu) \
         + U3/(2*pi) * chinnery(uxx_tf,x,p,L,W,q,delta,nu)
    uxy = -U1/(2*pi) * chinnery(uxy_ss,x,p,L,W,q,delta,nu) \
         - U2/(2*pi) * chinnery(uxy_ds,x,p,L,W,q,delta,nu) \
         + U3/(2*pi) * chinnery(uxy_tf,x,p,L,W,q,delta,nu)
    uyx = -U1/(2*pi) * chinnery(uyx_ss,x,p,L,W,q,delta,nu) \
         - U2/(2*pi) * chinnery(uyx_ds,x,p,L,W,q,delta,nu) \
         + U3/(2*pi) * chinnery(uyx_tf,x,p,L,W,q,delta,nu)
    uyy = -U1/(2*pi) * chinnery(uyy_ss,x,p,L,W,q,delta,nu) \
         - U2/(2*pi) * chinnery(uyy_ds,x,p,L,W,q,delta,nu) \
         + U3/(2*pi) * chinnery(uyy_tf,x,p,L,W,q,delta,nu)
    
    # Rotation from Okada's axes to geographic
    unn = sin(strike)**2*uxx + sin(2*strike)*(uxy + uyx)/2.0 + cos(strike)**2*uyy
    une = sin(2*strike)*(uxx - uyy)/2.0 + cos(strike)**2*uyx - sin(strike)**2*uxy
    uen = sin(2*strike)*(uxx - uyy)/2.0 - sin(strike)**2*uyx + cos(strike)**2*uxy
    uee = cos(strike)**2*uxx - sin(2*strike)*(uyx + uxy)/2.0 + sin(strike)**2*uyy
    
    # Assigns output arguments

    return (ue,un,uz,uze,uzn,unn,une,uen,uee)
    
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
# Notes for I... and K... subfunctions:
#
#    1. original formulas use Lame's constant as mu/(mu+lambda) which
#       depends only on the Poisson's ratio = 1 - 2*nu
#    2. tests for cos(delta) == 0 are made with "cos(delta) > eps" 
#       because cos(90*pi/180) is not zero but = 6.1232e-17 (!)
    
    
# =================================================================
# Chinnery's notation [equation (24) p. 1143]
    
# -----------------------------------------------------------------
def chinnery(f,x,p,L,W,q,delta,nu):
     return f(x,p,q,delta,nu) \
          - f(x,p-W,q,delta,nu) \
          - f(x-L,p,q,delta,nu) \
          + f(x-L,p-W,q,delta,nu)
    
    
# =================================================================
# Displacement subfunctions
    
# strike-slip displacement subfunctions [equation (25) p. 1144]
    
# -----------------------------------------------------------------
def ux_ss(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    # pdb.set_trace()
    return xi*q/(R*(R + eta)) \
        + atan((xi*eta)/(q*R)) \
        + I1(xi,eta,q,delta,nu,R)*sin(delta)
    
# -----------------------------------------------------------------
def uy_ss(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    # pdb.set_trace()
    return (eta*cos(delta) + q*sin(delta))*q/(R*(R + eta)) \
        + q*cos(delta)/(R + eta) \
        + I2(eta,q,delta,nu,R)*sin(delta)
    
# -----------------------------------------------------------------
def uz_ss(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(delta) - q*cos(delta)
    # pdb.set_trace()
    return (eta*sin(delta) - q*cos(delta))*q/(R*(R + eta)) \
        + q*sin(delta)/(R + eta) \
        + I4(db,eta,q,delta,nu,R)*sin(delta)
    
# dip-slip displacement subfunctions [equation (26) p. 1144]
# -----------------------------------------------------------------
def ux_ds(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return q/R \
        - I3(eta,q,delta,nu,R)*sin(delta)*cos(delta)
    
# -----------------------------------------------------------------
def uy_ds(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return (eta*cos(delta) + q*sin(delta))*q/(R*(R + xi)) \
        + cos(delta)*atan((xi*eta)/(q*R)) \
        - I1(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta)
    
# -----------------------------------------------------------------
def uz_ds(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(delta) - q*cos(delta)
    return db*q/(R*(R + xi)) \
        + sin(delta)*atan((xi*eta)/(q*R)) \
        - I5(xi,eta,q,delta,nu,R,db)*sin(delta)*cos(delta)
    
# tensile fault displacement subfunctions [equation (27) p. 1144]
# -----------------------------------------------------------------
def ux_tf(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return q**2 /(R*(R + eta)) \
        - I3(eta,q,delta,nu,R)*sin(delta)**2
    
# -----------------------------------------------------------------
def uy_tf(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return -(eta*sin(delta) - q*cos(delta))*q/(R*(R + xi)) \
        - sin(delta)*(xi*q/(R*(R + eta)) \
        - atan((xi*eta)/(q*R))) \
        - I1(xi,eta,q,delta,nu,R)*sin(delta)**2
    
# -----------------------------------------------------------------
def uz_tf(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(delta) - q*cos(delta)
    return (eta*cos(delta) + q*sin(delta))*q/(R*(R + xi)) \
        + cos(delta)*(xi*q/(R*(R + eta)) \
        - atan((xi*eta)/(q*R))) \
        - I5(xi,eta,q,delta,nu,R,db)*sin(delta)**2
    
    
# I... displacement subfunctions [equations (28) (29) p. 1144-1145]
# -----------------------------------------------------------------
def I1(xi,eta,q,delta,nu,R):
    db = eta*sin(delta) - q*cos(delta)
    if cos(delta) > eps:
        return (1 - 2*nu) * (-xi/(cos(delta)*(R+db))) - sin(delta)/cos(delta) *I5(xi,eta,q,delta,nu,R,db)
    else:
        return -(1 - 2*nu)/2 * xi*q/(R + db)**2
    
# -----------------------------------------------------------------
def I2(eta,q,delta,nu,R):
    return (1 - 2*nu) * (-log(R + eta)) - I3(eta,q,delta,nu,R)
    
# -----------------------------------------------------------------
def I3(eta,q,delta,nu,R):
    yb = eta*cos(delta) + q*sin(delta)
    db = eta*sin(delta) - q*cos(delta)
    if cos(delta) > eps:
        return (1 - 2*nu) * (yb/(cos(delta)*(R + db)) - log(R + eta)) \
            + sin(delta)/cos(delta) * I4(db,eta,q,delta,nu,R)
    else:
        return (1 - 2*nu)/2 * (eta/(R + db) + yb*q/(R + db)**2 - log(R + eta))
    
# -----------------------------------------------------------------
def I4(db,eta,q,delta,nu,R):
    if cos(delta) > eps:
        return (1 - 2*nu) * 1/cos(delta) * (log(R + db) - sin(delta)*log(R + eta))
    else:
        return -(1 - 2*nu) * q/(R + db)
    
# -----------------------------------------------------------------
def I5(xi,eta,q,delta,nu,R,db):
    X = sqrt(xi**2 + q**2)
    if cos(delta) > eps:
        return (1 - 2*nu) * 2/cos(delta) \
            * atan((eta*(X + q*cos(delta)) + X*(R + X)*sin(delta))/(xi*(R + X)*cos(delta)))
    else:
        return -(1 - 2*nu) * xi*sin(delta)/(R + db)
    
    
# =================================================================
# Tilt subfunctions
    
# strike-slip tilt subfunctions [equation (37) p. 1147]
    
# -----------------------------------------------------------------
def uzx_ss(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return -xi*q**2*A(eta,R)*cos(delta) \
        + ((xi*q)/R**3 - K1(xi,eta,q,delta,nu,R))*sin(delta)
    
# -----------------------------------------------------------------
def uzy_ss(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(delta) - q*cos(delta)
    yb = eta*cos(delta) + q*sin(delta)
    return (db*q/R**3)*cos(delta) \
        + (xi**2*q*A(eta,R)*cos(delta) - sin(delta)/R + yb*q/R**3 - K2(xi,eta,q,delta,nu,R))*sin(delta)
    
# dip-slip tilt subfunctions [equation (38) p. 1147]
    
# -----------------------------------------------------------------
def uzx_ds(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(delta) - q*cos(delta)
    return db*q/R**3 \
        + q*sin(delta)/(R*(R + eta)) \
        + K3(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta)
    
# -----------------------------------------------------------------
def uzy_ds(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(delta) - q*cos(delta)
    yb = eta*cos(delta) + q*sin(delta)
    return yb*db*q*A(xi,R) \
        - (2*db/(R*(R + xi)) + xi*sin(delta)/(R*(R + eta)))*sin(delta) \
        + K1(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta)
    
# tensile fault tilt subfunctions [equation (39) p. 1147]
    
# -----------------------------------------------------------------
def uzx_tf(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return q**2/R**3*sin(delta) \
        - q**3*A(eta,R)*cos(delta) \
        + K3(xi,eta,q,delta,nu,R)*sin(delta)**2
    
# -----------------------------------------------------------------
def uzy_tf(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(delta) - q*cos(delta)
    yb = eta*cos(delta) + q*sin(delta)
    return (yb*sin(delta) + db*cos(delta))*q**2*A(xi,R) \
        + xi*q**2*A(eta,R)*sin(delta)*cos(delta) \
        - (2*q/(R*(R + xi)) - K1(xi,eta,q,delta,nu,R))*sin(delta)**2
    
# -----------------------------------------------------------------
def A(x,R):
    return (2*R + x)/(R**3*(R + x)**2)
    
# K... tilt subfunctions [equations (40) (41) p. 1148]
# -----------------------------------------------------------------
def K1(xi,eta,q,delta,nu,R):
    db = eta*sin(delta) - q*cos(delta)
    if cos(delta) > eps:
        return (1 - 2*nu) * xi/cos(delta) * (1/(R*(R + db)) - sin(delta)/(R*(R + eta)))
    else:
        return (1 - 2*nu) * xi*q/(R + db)**2
    
# -----------------------------------------------------------------
def K2(xi,eta,q,delta,nu,R):
    return (1 - 2*nu) * (-sin(delta)/R + q*cos(delta)/(R*(R + eta))) - K3(xi,eta,q,delta,nu,R)
    
# -----------------------------------------------------------------
def K3(xi,eta,q,delta,nu,R):
    db = eta*sin(delta) - q*cos(delta)
    yb = eta*cos(delta) + q*sin(delta)
    if cos(delta) > eps:
        return (1 - 2*nu) * 1/cos(delta) * (q/(R*(R + eta)) - yb/(R*(R + db)))
    else:
        return (1 - 2*nu) * sin(delta)/(R + db) * (xi**2/(R*(R + db)) - 1)
    
    
# =================================================================
# Strain subfunctions
    
# strike-slip strain subfunctions [equation (31) p. 1145]
    
# -----------------------------------------------------------------
def uxx_ss(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return xi**2*q*A(eta,R) \
        - J1(xi,eta,q,delta,nu,R)*sin(delta)
    
# -----------------------------------------------------------------
def uxy_ss(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(delta) - q*cos(delta)
    return xi**3*db/(R**3*(eta**2 + q**2)) \
        - (xi**3*A(eta,R) + J2(xi,eta,q,delta,nu,R))*sin(delta)
    
# -----------------------------------------------------------------
def uyx_ss(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return xi*q/R**3*cos(delta) \
        + (xi*q**2*A(eta,R) - J2(xi,eta,q,delta,nu,R))*sin(delta)
    
# -----------------------------------------------------------------
def uyy_ss(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    yb = eta*cos(delta) + q*sin(delta)
    return yb*q/R**3*cos(delta) \
        + (q**3*A(eta,R)*sin(delta) - 2*q*sin(delta)/(R*(R + eta)) \
            - (xi**2 + eta**2)/R**3*cos(delta) - J4(xi,eta,q,delta,nu,R))*sin(delta)
        
# dip-slip strain subfunctions [equation (32) p. 1146]
    
# -----------------------------------------------------------------
def uxx_ds(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return xi*q/R**3 \
        + J3(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta)
    
# -----------------------------------------------------------------
def uxy_ds(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    yb = eta*cos(delta) + q*sin(delta)
    return yb*q/R**3 \
        - sin(delta)/R \
        + J1(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta)
    
# -----------------------------------------------------------------
def uyx_ds(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    yb = eta*cos(delta) + q*sin(delta)
    return yb*q/R**3 \
        + q*cos(delta)/(R*(R + eta)) \
        + J1(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta)
    
# -----------------------------------------------------------------
def uyy_ds(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    yb = eta*cos(delta) + q*sin(delta)
    return yb**2*q*A(xi,R) \
        - (2*yb/(R*(R + xi)) + xi*cos(delta)/(R*(R + eta)))*sin(delta) \
        + J2(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta)
    
# tensile fault strain subfunctions [equation (33) p. 1146]
    
# -----------------------------------------------------------------
def uxx_tf(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return xi*q**2*A(eta,R) \
        + J3(xi,eta,q,delta,nu,R)*sin(delta)**2
    
# -----------------------------------------------------------------
def uxy_tf(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(delta) - q*cos(delta)
    return -db*q/R**3 \
        - xi**2*q*A(eta,R)*sin(delta) \
        + J1(xi,eta,q,delta,nu,R)*sin(delta)**2
    
# -----------------------------------------------------------------
def uyx_tf(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    return q**2/R**3*cos(delta) \
        + q**3*A(eta,R)*sin(delta) \
        + J1(xi,eta,q,delta,nu,R)*sin(delta)**2
    
# -----------------------------------------------------------------
def uyy_tf(xi,eta,q,delta,nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(delta) - q*cos(delta)
    yb = eta*cos(delta) + q*sin(delta)
    return (yb*cos(delta) - db*sin(delta))*q**2*A(xi,R) \
        - q*sin(2*delta)/(R*(R + xi)) \
        - (xi*q**2*A(eta,R) - J2(xi,eta,q,delta,nu,R))*sin(delta)**2
    
    
# J... tensile fault subfunctions [equations (34) (35) p. 1146-1147]
# -----------------------------------------------------------------
def J1(xi,eta,q,delta,nu,R):
    db = eta*sin(delta) - q*cos(delta)
    if cos(delta) > eps:
        return (1 - 2*nu) * 1/cos(delta) * (xi**2/(R*(R + db)**2) - 1/(R + db)) \
            - sin(delta)/cos(delta)*K3(xi,eta,q,delta,nu,R)
    else:
        return (1 - 2*nu)/2 * q/(R + db)**2 * (2*xi**2/(R*(R + db)) - 1)
    
# -----------------------------------------------------------------
def J2(xi,eta,q,delta,nu,R):
    db = eta*sin(delta) - q*cos(delta)
    yb = eta*cos(delta) + q*sin(delta)
    if cos(delta) > eps:
        return (1 - 2*nu) * 1/cos(delta) * xi*yb/(R*(R + db)**2) \
            - sin(delta)/cos(delta)*K1(xi,eta,q,delta,nu,R)
    else:
        return (1 - 2*nu)/2 * xi*sin(delta)/(R + db)**2 * (2*q**2/(R*(R + db)) - 1)
    
# -----------------------------------------------------------------
def J3(xi,eta,q,delta,nu,R):
    return (1 - 2*nu) * -xi/(R*(R + eta)) \
        - J2(xi,eta,q,delta,nu,R)
    
# -----------------------------------------------------------------
def J4(xi,eta,q,delta,nu,R):
    return (1 - 2*nu) * (-cos(delta)/R - q*sin(delta)/(R*(R + eta))) \
        - J1(xi,eta,q,delta,nu,R)



if __name__ == "__main__":
    import sys
    import re
    if len(sys.argv) < 3:
        print "Syntax: input_file output_file"
        sys.exit()

    with open(sys.argv[1]) as infile:
        with open(sys.argv[2],"w") as outfile:
            outfile.write("x\ty\tux\tuy\tuz\n")
            line = infile.readline()
            (E0,N0,DEPTH,STRIKE,DIP,LENGTH,WIDTH,U1,U2,U3) = map(
                lambda(x): float(x), line.split())
            for l in infile:
                points = None
                try:
                    parts=map(lambda(x):float(x),l.split())
                    if len(parts) == 6:
                        nx = (int)(parts[4])
                        ny = (int)(parts[5])
                        points = [
                            ((parts[0]*(nx-i) + parts[2]*i)/nx,
                            (parts[1]*(ny-j) + parts[3]*j)/ny)
                            for i in range(nx+1)
                            for j in range(ny+1)
                            ]
                    else:
                        points = [map(lambda(x): float(x),l.split())]

                    for p in points: 
                        try:
                            (ue,un,uz,uze,uzn,unn,une,uen,uee) = okada85(
                               p[0]-E0,p[1]-N0,DEPTH,STRIKE,DIP,LENGTH,WIDTH,U1,U2,U3)
                            outfile.write("%.2f\t%.2f\t%.5f\t%.5f\t%.5f\n"%(p[0],p[1],ue,un,uz))
                        except:
                            sys.stderr.write("Error: " + str(sys.exc_info()[1]) + "\n")
                except:
                    sys.stderr.write("Error: " + str(sys.exc_info()[1]) + "\n")
