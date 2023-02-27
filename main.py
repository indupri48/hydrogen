from rtlsdr import RtlSdr
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

# Ephemeris
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, ICRS, EarthLocation
from datetime import datetime

def equatorial(lat, long, alt, az):

	TIME = Time(datetime.utcnow())
	QTH = EarthLocation(lat=lat*u.degree, lon=long*u.degree, height=0)

	horizontal_coord = AltAz(alt=alt*u.degree, az = az*u.degree, pressure = 0*u.bar, obstime=TIME, location=QTH)
	equatorial_coord = SkyCoord(horizontal_coord.transform_to(ICRS()))

	return round(equatorial_coord.ra.degree, 2), round(equatorial_coord.dec.degree, 2)

def measure_spectrum(sdr, N_SAMPLES, N_FFT):

	total_psd = np.zeros(N_SAMPLES)

	for i in range(N_FFT):

		samples = sdr.read_samples(N_SAMPLES)
		psd = np.abs(np.fft.fft(samples)) ** 2 / (N_SAMPLES * fs)
		psd = 10.0*np.log10(psd)
		psd = np.fft.fftshift(psd)

		total_psd = np.add(total_psd, psd)

	mean_psd = np.true_divide(total_psd, N_FFT)

	return mean_psd

def take_reading(sdr, long, lat, alt, az):

	ra, dec = equatorial(lat, long, alt, az)

	psd = measure_spectrum(sdr, 1024, 10000)

	frequency = np.arange(-sdr.sample_rate/2.0, sdr.sample_rate/2.0, sdr.sample_rate/len(psd))
	frequency += sdr.center_freq

	peak_frequency = frequency[np.argmax(psd)]

	print(ra, dec, peak_frequency)

	'''
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Magnitude (dB)')
	plt.grid(True)
	plt.plot(frequency, psd)
	plt.axvline(peak_frequency, color='black', linestyle='dashed')
	plt.show()
	'''
fs = 2.048e6
cf = 101e6

# Initialize SDR object
sdr = RtlSdr()
sdr.center_freq = cf
sdr.sample_rate = fs
sdr.freq_correction = 60
sdr.gain = 'auto'

lat = input('Enter your latitude (-180 to 180 degrees): ')
long = input('Enter your longitude (-90 to 90 degrees): ')

while True:

	alt = input('Enter altitude (0 to 90 degrees): ')
	az = input('Enter azimuth (0 to 360 degrees): ')

	take_reading(sdr, lat, long, alt, az)


'''
plt.xlabel('Frequency (MHz)')
plt.ylabel('Magnitude (dB)')
plt.grid(True)
plt.plot(frequency, psd)
plt.axvline(peak_frequency, color='black', linestyle='dashed')
plt.show()
'''
