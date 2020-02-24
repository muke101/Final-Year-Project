import matplotlib.pyplot as plt
import csv
import sys

results = csv.reader(open(str(sys.argv[1]), 'r'))
x1 = []
y1 = []
x2 = []
y2 = []

def format(res):
    return float(res.split(':')[1].strip())

for row in results:
    x1.append(format(row[1]))
    y1.append(format(row[2]))
    if row[0] == '20':
        break

for row in results:
    x2.append(format(row[1]))
    y2.append(format(row[2]))

fig = plt.figure()
ax = fig.add_subplot()
ax.plot(x1,y1,label="<ln(fc)>nf")
ax.plot(x2,y2,label="<ln(fc)>Ca")
ax.set_xlabel('x')
ax.legend()
plt.show()
