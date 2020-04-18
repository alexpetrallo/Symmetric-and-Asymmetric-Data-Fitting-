import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

#import iodine concentration data
data = np.loadtxt('iodine.txt')
time = data[:,0]
concentration = data[:,1]

#plotting time vs concentration
plt.scatter(time, concentration)
plt.xlabel('Time (s)')
plt.ylabel('Iodine Concentration (Geiger Counts/s)')
plt.title('Iodine Concentration in Blood over Time')
plt.savefig('timeVconc')
# plt.show()

#taking natural log of time and concentration
ln_time = np.log(time)
ln_conc = np.log(concentration)
'''
#plotting natural log of time and natural log of concentration
fig2 = plt.figure(2)
plt.scatter(ln_time, ln_conc)
plt.xlabel('Natural Log of Time (s)')
plt.ylabel('Natural Log of Iodine Concentration (Geiger Counts/s)')
plt.title('Ln time vs Ln concentration')
plt.show()
'''
#linearizing
xn = (4.5, 8, 200)
m, c = np.polyfit(ln_time, ln_conc, 1)
yn = np.polyval([m, c], xn)

#plotting natural log of time & natural log of concentration w/ fit
fig2 = plt.figure(2)
plt.scatter(ln_time, ln_conc)
plt.plot(xn, yn, 'k')
plt.xlim(4, 8)
plt.ylim(7, 9.5)
plt.xlabel('Natural Log of Time (s)')
plt.ylabel('Natural Log of Iodine Concentration (Geiger Counts/s)')
plt.title('Ln time vs Ln concentration with fit')
plt.savefig('lnTimeVconc')
# plt.show()

#linear substitution
a = np.exp(c)
b = m
def exp(x):
    a = np.exp(c)
    b = m
    return a*np.exp(b*x)

#calculate exponential fit & r2
expfit = exp(ln_time)
r2 = r2_score(concentration, expfit)
print ('The r2 for the simple exponential fit is:')
print(r2)

#plotting time and concentration w/ fit
fig3 = plt.figure(3)
plt.plot(time, expfit, 'k')
plt.scatter(time, concentration)
plt.xlabel('Time (s)')
plt.ylabel('Iodine Concentration (Geiger Counts/s)')
plt.title('Time vs Concentration with fit. $R^2$ =%.5f' % r2)
plt.savefig('tVcWfit')
plt.show()

midRangeTime = [40,]
midRangeConc = [12089.5,]
for i in range(len(time)-1):
    midRangeTime.append((time[i] + time[i+1])/2)
for i in range(len(concentration)-1):
    midRangeConc.append((concentration[i] + concentration[i+1])/2)

ln_time = np.log(time)
ln_conc = np.log(concentration)

midtime = midRangeTime
midconc = midRangeConc



#square-root times
sqrt_time = np.sqrt(time)
sqrt_midtime = np.sqrt(midtime)

#linearizing sqrt times
xn3 = (10, 50, 20)
m3, c3 = np.polyfit(sqrt_time, ln_conc, 1)
yn3 = np.polyval([m3, c3], xn3)

#linearizing sqrt mid-times
m3m, c3m = np.polyfit(sqrt_midtime, ln_conc, 1)
yn3m = np.polyval([m3m, c3m], xn3)

#plotting sqrt time w/ natural log concentration
fig4 = plt.figure(4)
plt.scatter(sqrt_time, ln_conc)
plt.plot(xn3, yn3, 'k')
plt.xlim(0, 50)
plt.ylim(7, 10)
plt.xlabel('$\sqrt{t}$ (s)')
plt.ylabel('Iodine Ln Concentration (Geiger Counts/s)')
plt.title('$\sqrt{Time}$ vs Ln Concentration with linear fit')
plt.savefig('sqrtTime')
plt.show()

#linear susbstitution w/ sqrt times
a3 = np.exp(c3)
b3 = m3
def exp3(x):
    a3 = np.exp(c3)
    b3 = m3
    return a3*np.exp(b3*np.sqrt(x))

#calculate exponential fit & r2
expfit3 = exp3(time)
r23 = r2_score(concentration, expfit3)
print ('The r2 for the sqrt t exponential fit is:')
print(r23)

#plotting sqrt time and concentration
fig5 = plt.figure(5)
plt.plot(time, expfit3, 'k')
plt.scatter(time, concentration)
plt.xlabel('Time (s)')
plt.ylabel('Iodine Concentration (Geiger Counts/s)')
plt.title('Time vs Concentration with $e^{- \lambda \sqrt{t}}$ fit. $R^2$ =%.5f' % r23)
plt.savefig('lambdaFit')
plt.show()

#plotting sqrt mid-times and natural log concentrations
fig6 = plt.figure(6)
plt.scatter(sqrt_midtime, ln_conc)
plt.plot(xn3, yn3m, 'k')
plt.xlim(0, 50)
plt.ylim(7, 10)
plt.xlabel('Time (s)')
plt.ylabel('Iodine Concentration (Geiger Counts/s)')
plt.title('Mid-range time vs Ln concentration with linear fit')
plt.savefig('midRangeLinFit')
plt.show()

#linear susbstitution w/ sqrt mid-times
a3m = np.exp(c3m)
b3m = m3m
def exp3m(x):
    a3m = np.exp(c3m)
    b3m = m3m
    return a3m*np.exp(b3m*np.sqrt(x))

#calculate exponential fit & r2
expfit3m = exp3m(midtime)
r23m = r2_score(concentration, expfit3m)
print ('The r2 for the sqrt t exponential fit w/ mid-range times is:')
print(r23m)

#plotting midtimes & concentrations w/ fit
fig7 = plt.figure(7)
plt.plot(midtime, expfit3m, 'k')
plt.scatter(midtime, concentration)
plt.xlabel('Time (s)')
plt.ylabel('Iodine Concentration (Geiger Counts/s)')
plt.title('Mid-range Time vs Concentration with $e^{- \lambda \sqrt{t}}$ fit. $R^2$ =%.5f' % r23m)
plt.show()
plt.savefig('midRangeExpFit')

e3 = exp3m(midtime)
r23mW = r2_score(midconc, e3)
fig8 = plt.figure(8)
plt.plot(midtime, e3, 'k')
plt.scatter(midtime, midconc)
plt.xlabel('Time (s)')
plt.ylabel('Iodine Concentration (Geiger Counts/s)')
plt.title('Mid-range Time vs Mid-range Concentration with $e^{- \lambda \sqrt{t}}$ fit. $R^2$ =%.5f' % r23mW)
plt.savefig('midRangeExpFitWorse')


constants = np.linspace(0, 9999, 10000)


expfit_const = []
r2_list = []
r2_best = 0
for i in constants:
    exp_new = exp3(time) + i
    expfit_const = np.append(expfit_const, exp_new)
    r2_list = np.append(r2_list, r2_score(concentration, exp_new))
    r2_const = r2_score(concentration, exp_new)
    if r2_const > r2_best:
        r2_best = r2_const
print(r2_best)
# print('The additive constant is', constants[index])
index = np.argmax(r2_list)

fig9 = plt.figure(9)
plt.scatter(time, concentration)
plt.plot(time, exp3(time) + constants[index], 'k')
plt.title('Iodine Purification Data fit with $e^{- \lambda \sqrt{t}} + c$ $R^2$ =%.5f' % r2_best)
plt.savefig('expAddConst')
# plt.show()
