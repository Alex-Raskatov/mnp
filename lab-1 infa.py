import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# The method of least squares for the linear dependence
def MLS(x, y):
    l = len(x)

    a = (sum(x * y) / l - sum(x) * sum(y) / (l ** 2)) / (sum(x ** 2) / l - sum(x) ** 2 / (l ** 2))
    b = sum(y) / l - a * sum(x) / l

    s_a = (1 / l ** 0.5) * ((sum(y ** 2) / l - sum(y) ** 2 / l ** 2) / (sum(x ** 2) / l - sum(x) ** 2 / l ** 2) - a ** 2) ** 0.5
    s_b = s_a * (sum(x ** 2) / l - sum(x) ** 2 / l ** 2) ** 0.5
    return a, b, s_a, s_b

# Data reading from input.txt (data must be divided by spacebar).
# The first string must be x, the second must be y.
file = open('./data/task1/out_000.dat', 'r', encoding='utf-8')
t = list(map(lambda x: x.rstrip().split(), file.readlines()))

x = []
y = []

for i in range(len(t)):
    x.append(i)
    y.append(float(t[i][0]))

x = np.array(x)
y = np.array(y)

#x = np.log(x)
#y = np.log(y)

if len(x) != len(y):
    print('ERROR: len(x) != len(y)')
else:
    # Setting the parameters of the graph
    plt.figure(figsize=(10,6))

    # Setting the limits of axes (you should set it later be yourself)
    x_min = min(x)
    x_max = max(x)
    y_min = min(y)
    y_max = max(y)
    plt.ylim([y_min - (y_max-y_min)*0.05, y_max + (y_max-y_min)*0.05])
    plt.xlim([x_min - (x_max-x_min)*0.05, x_max + (x_max-x_min)*0.05])

    # Plotting the graphs
    a, b, s_a, s_b = MLS(x, y)
    x_MLS = np.linspace(x_min - 0.05*(x_max-x_min), x_max + 0.05*(x_max-x_min), 2)
    y_MLS = a*x_MLS+b
    plt.plot(x_MLS, y_MLS, color = 'grey')
    plt.plot(x, y,'o', color = 'black')
    plt.title('')
    plt.ylabel('')
    plt.xlabel('')
    plt.minorticks_on()
    plt.grid(visible=True, which='major', axis='both', alpha=1)
    plt.grid(visible=True, which='minor', axis='both', alpha=0.2)
    print('y = ax + b \n')
    print(f'a = {a}')
    print(f'b = {b}')
    print(f'sigma_a = {s_a}')
    print(f'sigma_b = {s_b}')
    print(f'epsilon_a = {abs(int(s_a/a * 1000)/10)}%')
    print(f'epsilon_b = {abs(int(s_b/b * 1000)/10)}% \n')
    print(f'x_min = {x_min}')
    print(f'x_max = {x_max}')
    print(f'y_min = {y_min}')
    print(f'y_max = {y_max} \n')
    plt.show()
    print('Thank you for using the app!')