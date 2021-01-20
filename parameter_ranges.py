import numpy as np

directory = "/srv/djw1g16/cloudy_tables/test_tables"

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

dir_format = '{directory}/n{n:.0f}/In{i:.1f}/Te{t:.2f}/'
file_format = "n{n:.0f}_In{i:.1f}_Te{t:.2f}"

def format_dir(n,i,t):
    return dir_format.format(directory=directory,n=n,i=i,t=t)

def format_file(n,i,t):
    return file_format.format(n=n,i=i,t=t)
