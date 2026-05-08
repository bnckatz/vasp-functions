# There are obviously lots of PROCAR parsers out there, but I eventually found myself wanting my own that read the data into a
# simply interpretable structure.  I vibecoded the basic structure, then tore the guts out and rewrote most of it.
# My thanks to whoever's code I'm sure the AI ingested to create the basic structure in the first place...

import numpy as np
import re

def read_procar(filename):
    """
    Minimal PROCAR parser for spin-polarized (no SOC) VASP runs.
    Returns structured numpy arrays suitable for Fermi surface plotting.
    """

    with open(filename, 'r') as f:
        lines = f.readlines()

    # --- Parse header ---
    line_index = 1
    header_match = re.search(
        r'# of k-points:\s+(\d+)\s+# of bands:\s+(\d+)\s+# of ions:\s+(\d+)',
        lines[line_index])

    if header_match is None:
        raise ValueError("Could not parse PROCAR header.")

    n_kpoints = int(header_match.group(1))
    n_bands = int(header_match.group(2))
    n_ions = int(header_match.group(3))

    # Detect spin-polarization
    count = 0
    for line in lines:
        if re.findall(r'ions',line)!=[]:
            count += 1
    if count ==2:
        print("Spin-polarization present")
        n_spin = 2
    else:
        print("Spin-polarization absent")
        n_spin = 1

    kpoints = np.zeros((n_kpoints, 3)) # x,y,z kpoint coords
    weights = np.zeros(n_kpoints) # a weight for each kpoint
    energies = np.zeros((n_spin, n_kpoints, n_bands)) # energies for state sorted by [spin, kpoint #, band #]
    occupations = np.zeros((n_spin, n_kpoints, n_bands)) # occupations for state sorted by [spin, kpoint #, band #]
    # find which projections are included:
    line_index = 7
    text_match=re.findall(r'\w*-?\w+',lines[line_index])
    text_match=text_match[1:len(text_match)] # strip first entry, this gives us our number of projections (s,p...tot)
    # check if the PROCAR is noncollinear or not
    wordfind=re.compile(r'\w+') # compile this for re-use
    while (wordfind.findall(lines[line_index])==[]
        or wordfind.findall(lines[line_index])[0]!="tot"):
        line_index += 1
    line_index += 1
    if wordfind.findall(lines[line_index])==[]:
        print("PROCAR spin is collinear")
        print(lines[line_index])
        print(wordfind.findall(lines[line_index]))
        ncollin=False
        projections=np.zeros((n_ions,n_kpoints,n_bands,len(text_match))) # projections sorted by ion,kpoint,band
    else:
        print("PROCAR spin is non-collinear")
        ncollin=True
        projections = np.zeros((n_ions, n_kpoints, n_bands, len(text_match),3))


    # reset line index to begin processing
    line_index = 2
    # compile regex we're going to use a lot
    floatfind = re.compile(r'\d+\.\d+')

    for spin in range(n_spin):
        print("Spin",spin)
        for k in range(n_kpoints):
            # Skip empty lines
            while lines[line_index].strip() == "":
                line_index += 1

            # k-point line
            # Example:
            # k-point   1 :    0.00000000 0.00000000 0.00000000     weight = 0.12500000
            if spin==1:  # only need to assign on first pass through
                kp_match = floatfind.findall(lines[line_index])
                kpoints[k] = [
                    float(kp_match[0]),
                    float(kp_match[1]),
                    float(kp_match[2])
                ]
                weights[k] = float(kp_match[3])
            line_index += 1

            for band in range(n_bands):
                # find first band line
                while (re.findall(r'\w+',lines[line_index])==[] # if it's empty
                       or re.findall(r'\w+',lines[line_index])[0]!="band"): # if it's the 'total' proj line
                    line_index += 1

                # band line
                # Example:
                # band  1 # energy  -5.12345678 # occ.  1.00000000
                band_match = floatfind.findall(lines[line_index])

                energies[spin, k, band] = float(band_match[0])
                occupations[spin, k, band] = float(band_match[1])
                line_index += 1

                # Skip to projection block:
                line_index += 1
                for i in range(1+3*ncollin):
                    line_index += 1
                    for ion in range(n_ions):
                        # match projections
                        if ncollin==True:
                            projections[ion,k,band,:,int(np.floor(ion/n_ions))] = floatfind.findall(lines[line_index])
                        else:
                            projections[ion,k,band,:]=floatfind.findall(lines[line_index])
                        line_index += 1

    return {
        "n_kpoints": n_kpoints,
        "n_bands": n_bands,
        "n_ions": n_ions,
        "kpoints": kpoints,
        "weights": weights,
        "energies": energies,
        "occupations": occupations,
        "projections": projections,
        "proj_included": text_match
    }
    
    
