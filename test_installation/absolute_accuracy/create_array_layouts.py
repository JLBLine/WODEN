import numpy as np

def calc_b(phi_simple, num_mult):
    """Given the target angle `phi_simple` (radians), and
    an integer multiplier `num_mult`, find a baseline
    length `b` that results in the same sine/cosine
    angle"""

    if phi_simple == 0:

        b = 2*np.pi*num_mult

    else:
        b = (phi_simple + 2*np.pi*num_mult) / phi_simple

    return b

##Known list of angles that have predictable sin/cos outputs
phi_simples = [0.0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, 2*np.pi/3,
               3*np.pi/4, 5*np.pi/6, np.pi, 7*np.pi/6, 5*np.pi/4]

##For known angles, calculate l,m,n coords and check they are legit
for ind, phi_simple in enumerate(phi_simples):

    for num_mult in [1, 10, 100, 1000, 10000]:

        b = calc_b(phi_simple, num_mult)

        with open(f"array_layouts/array_layout_ang{ind:02d}_n{num_mult:05d}.txt", 'w') as outfile:
            outfile.write(f"{b} {b} {b}\n")
            outfile.write("0.0 0.0 0.0")

        # phi_outcome = 2*np.pi*b*(l + m + n - 1)
