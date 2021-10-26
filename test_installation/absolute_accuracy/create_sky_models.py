import numpy as np

def find_lmn(phi_simple):
    """Given phi_simple (radians), find a set of l,m,n coords
    that produce phi_simple in the measurement equation"""
    numer = np.sqrt(2)*np.sqrt(-phi_simple**2 + -4*np.pi*phi_simple + 8*np.pi**2) + phi_simple + 2*np.pi
    denom = 6*np.pi

    n = numer / denom

    l = np.sqrt(1 - n*n) / np.sqrt(2)
    m = l

    return l,m,n

##Known list of angles that have predictable sin/cos outputs
phi_simples = [0.0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, 2*np.pi/3,
               3*np.pi/4, 5*np.pi/6, np.pi, 7*np.pi/6, 5*np.pi/4]

all_ras = []
all_decs = []

##For known angles, calculate l,m,n coords and check they are legit
for ind, phi_simple in enumerate(phi_simples):

    l,m,n = find_lmn(phi_simple)

    angular_offset = np.arcsin(l)

    ##RA is stored in hours (why did I do this)

    dec = angular_offset * (180.0/np.pi)

    ## This maths should produce a ra,dec that should make the correct
    ## l,m for a latitude of zero
    ra = np.arcsin((np.sin(angular_offset) / np.cos(angular_offset))) * (12.0/np.pi)

    all_ras.append(ra)
    all_decs.append(dec)

    with open(f"sky_models/srclist_ang{ind:02d}.txt", 'w') as outfile:
        outfile.write("SOURCE testing P 1 G 0 S 0 0\n")
        outfile.write(f"COMPONENT POINT {ra:.16f} {dec:.16f}\n")
        outfile.write(f"LINEAR 150e+6 1.0 0 0 0 0.0\n")
        outfile.write(f"ENDCOMPONENT\n")
        outfile.write(f"ENDSOURCE")
