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
