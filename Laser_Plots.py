from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

rcParams["font.family"] = "sans-serif"
rcParams['font.sans-serif'] = ['Tahoma']


def voltToWatt(volt, readout):
    return (volt/500)*readout

def average(lower, upper):
    return (lower + upper)/2
    
def error(lower, upper):
    return (upper - lower)
    
def ratioTransmitted(before, trans):
    return (trans/before)
    
def multError(laserAvg, transmittedAvg, laserErr, transmittedErr, R):
    return (R * np.sqrt(np.power((laserErr/laserAvg),2)+np.power((transmittedErr/transmittedAvg),2)))                    

distancetest_input = np.loadtxt("distancetest.txt")
distancetest = np.transpose(distancetest_input)

distancetest_dist = distancetest[0]
distancetest_min_mv = distancetest[1]
distancetest_max_mv = distancetest[2]

distancetest_min_mw = voltToWatt(distancetest_min_mv, 20)
distancetest_max_mw = voltToWatt(distancetest_max_mv, 20)

distancetest_avg_mw = average(distancetest_min_mw, distancetest_max_mw)
distancetest_err_mw = error(distancetest_min_mw, distancetest_max_mw)


conc100input = np.loadtxt("conc100.txt")
conc100 = np.transpose(conc100input)

conc100_length = conc100[0]
conc100_min_mv = conc100[1]
conc100_max_mv = conc100[2]
conc100_laser_avg = 13.448
conc100_laser_err = 0.064

conc100_min_mw = voltToWatt(conc100_min_mv, 20)
conc100_max_mw = voltToWatt(conc100_max_mv, 20)

conc100_avg_mw = average(conc100_min_mw, conc100_max_mw)
conc100_err_mw = error(conc100_min_mw, conc100_max_mw)

conc100_ratio_transmitted = ratioTransmitted(conc100_laser_avg, conc100_avg_mw)
conc100_ratio_transmitted_err = multError(conc100_laser_avg, conc100_avg_mw, conc100_laser_err, conc100_err_mw, conc100_ratio_transmitted)


conc050input = np.loadtxt("conc050.txt")
conc050 = np.transpose(conc050input)

conc050_length = conc050[0]
conc050_min_mv = conc050[1]
conc050_max_mv = conc050[2]
conc050_laser_avg = 13.59
conc050_laser_err = 0.26

conc050_min_mw = voltToWatt(conc050_min_mv, 20)
conc050_max_mw = voltToWatt(conc050_max_mv, 20)

conc050_avg_mw = average(conc050_min_mw, conc050_max_mw)
conc050_err_mw = error(conc050_min_mw, conc050_max_mw)

conc050_ratio_transmitted = ratioTransmitted(conc050_laser_avg, conc050_avg_mw)
conc050_ratio_transmitted_err = multError(conc050_laser_avg, conc050_avg_mw, conc050_laser_err, conc050_err_mw, conc050_ratio_transmitted)


conc025input = np.loadtxt("conc025.txt")
conc025 = np.transpose(conc025input)

conc025_length = conc025[0]
conc025_min_mv = conc025[1]
conc025_max_mv = conc025[2]
conc025_laser_avg = 13.564
conc025_laser_err = 0.336

conc025_min_mw = voltToWatt(conc025_min_mv, 20)
conc025_max_mw = voltToWatt(conc025_max_mv, 20)

conc025_avg_mw = average(conc025_min_mw, conc025_max_mw)
conc025_err_mw = error(conc025_min_mw, conc025_max_mw)

conc025_ratio_transmitted = ratioTransmitted(conc025_laser_avg, conc025_avg_mw)
conc025_ratio_transmitted_err = multError(conc025_laser_avg, conc025_avg_mw, conc025_laser_err, conc025_err_mw, conc025_ratio_transmitted)


conc000input = np.loadtxt("conc000.txt")
conc000 = np.transpose(conc000input)

conc000_length = conc000[0]
conc000_min_mv = conc000[1]
conc000_max_mv = conc000[2]
conc000_laser_avg = 13.606
conc000_laser_err = 0.316

conc000_min_mw = voltToWatt(conc000_min_mv, 20)
conc000_max_mw = voltToWatt(conc000_max_mv, 20)

conc000_avg_mw = average(conc000_min_mw, conc000_max_mw)
conc000_err_mw = error(conc000_min_mw, conc000_max_mw)

conc000_ratio_transmitted = ratioTransmitted(conc000_laser_avg, conc000_avg_mw)
conc000_ratio_transmitted_err = multError(conc000_laser_avg, conc000_avg_mw, conc000_laser_err, conc000_err_mw, conc000_ratio_transmitted)


dconcinput = np.loadtxt("dconc.txt")
dconc = np.transpose(dconcinput)

dconc_length = dconc[0]
dconc_min_mv = dconc[1]
dconc_max_mv = dconc[2]
dconc_laser_avg = 13.564
dconc_laser_err = 0.336

dconc_min_mw = voltToWatt(dconc_min_mv, 20)
dconc_max_mw = voltToWatt(dconc_max_mv, 20)

dconc_avg_mw = average(dconc_min_mw, dconc_max_mw)
dconc_err_mw = error(dconc_min_mw, dconc_max_mw)

dconc_ratio_transmitted = ratioTransmitted(dconc_laser_avg, dconc_avg_mw)
dconc_ratio_transmitted_err = multError(dconc_laser_avg, dconc_avg_mw, dconc_laser_err, dconc_err_mw, dconc_ratio_transmitted)


conc100_slope, conc100_intercept, conc100_r_value, conc100_p_value, conc100_std_err = stats.linregress(conc100_length, np.log(conc100_avg_mw))
conc050_slope, conc050_intercept, conc050_r_value, conc050_p_value, conc050_std_err = stats.linregress(conc050_length, np.log(conc050_avg_mw))
conc025_slope, conc025_intercept, conc025_r_value, conc025_p_value, conc025_std_err = stats.linregress(conc025_length, np.log(conc025_avg_mw))
dconc_slope, dconc_intercept, dconc_r_value, dconc_p_value, dconc_std_err = stats.linregress(dconc_length, np.log(dconc_avg_mw))


plt.figure(1)
plt.errorbar(distancetest_dist, distancetest_avg_mw, distancetest_err_mw, marker="x", linestyle="")
plt.tick_params(axis="both", direction="in")
plt.xlabel("Distance from Laser Aperature [cm]")
plt.ylabel("Power Transmitted [mW]")

plt.figure(2)
plt.errorbar(conc100_length, conc100_avg_mw, conc100_err_mw, marker="x", linestyle="", label="Concentration 1.00")
plt.errorbar(conc050_length*0.5, conc050_avg_mw, conc050_err_mw, marker="x", linestyle="", label="Concentration 0.50")
plt.errorbar(conc025_length*0.25, conc025_avg_mw, conc025_err_mw, marker="x", linestyle="", label="Concentration 0.25")

plt.tick_params(axis='both', direction="in")
plt.xlabel("Cell Length [cm]")
plt.ylabel("Power Transmitted [mW]")
plt.legend(loc='lower left')

plt.figure(3)
plt.errorbar(conc100_length, conc100_ratio_transmitted, conc100_ratio_transmitted_err, marker="x", linestyle="", label="Concentration 1.00")
plt.errorbar(conc050_length, conc050_ratio_transmitted, conc050_ratio_transmitted_err, marker=".", linestyle="", label="Concentration 0.50")
plt.errorbar(conc025_length, conc025_ratio_transmitted, conc025_ratio_transmitted_err, marker="P", linestyle="", label="Concentration 0.25")
plt.errorbar(conc000_length, conc000_ratio_transmitted, conc000_ratio_transmitted_err, marker="v", linestyle="", label="Concentration 0.00")

plt.tick_params(axis='both', direction="in")
plt.xlabel("Cell Length [cm]")
plt.ylabel("Ratio of Power Transmitted")
plt.legend(loc='lower left')

plt.figure(4)
plt.errorbar(dconc_length, dconc_ratio_transmitted, dconc_ratio_transmitted_err, marker="x", linestyle="")

plt.tick_params(axis='both', direction="in")
plt.xlabel("Cell Length * Concentration [number of particles]")
plt.ylabel("Ratio of Power Transmitted")
plt.semilogy()


plt.show()

