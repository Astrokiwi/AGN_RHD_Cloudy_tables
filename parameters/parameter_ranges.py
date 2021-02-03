import numpy as np

directory = "/srv/djw1g16/cloudy_tables/test_tables"

# radii of grain species 01-10, for both graphites and silicates (units=mm, I think?)
d_rads = np.array([6.080e-07,8.991e-07,1.330e-06,1.966e-06,2.907e-06,4.299e-06,6.358e-06,9.402e-06,1.390e-05,2.056e-05])

# Intensity at the sublimation radius
I0 = 5.6e7

# Intensities inside & outside rsub
logi = np.log10(I0) + np.linspace(2.5, -6.5, 19)

# Array of temperatures
logt = np.append(
    np.linspace(1, 5, 17)
    , np.linspace(6, 8, 3)
    )

logn = np.linspace(0, 7, 8)

table_tau = np.logspace(-2,np.log10(7),50)

dir_format = '{directory}/n{n:.0f}/In{i:.1f}/Te{t:.2f}/'
file_format = "n{n:.0f}_In{i:.1f}_Te{t:.2f}"

# key_order = ["n","i","t","arad","Htot","Cool","kabs","kscat","Tdust","dg"]
key_order = ["arad","Htot","Cool","kabs","kscat","Tdust","dg"]

def format_dir(n,i,t):
    return dir_format.format(directory=directory,n=n,i=i,t=t)

def format_file(*args):
    if len(args)==1:
        return format_file_rd(*args)
    else:
        return format_file_nit(*args)

def format_file_nit(n,i,t):
    return file_format.format(n=n,i=i,t=t)

def format_file_rd(rd):
    return file_format.format(n=rd.n,i=rd.i,t=rd.t)