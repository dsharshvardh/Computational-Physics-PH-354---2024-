# Question (a) - Plotting the given function 
import matplotlib.pyplot as plt

def P(x):
    return 924*(x**6) - 2772*(x**5) + 3150*(x**4) - 1680*(x**3) + 420*(x**2) - 42*x + 1

def deriv_P(x): # derivative of P(x)
    return 5544*(x**5) - 13860*(x**4) + 12600*(x**3) - 5040*(x**2) + 840*x -42

x_lst = [i*0.01 for i in range(0,101)]


P_lst = [P(x) for x in x_lst]

# Plotting
plt.plot(x_lst, P_lst, marker='o', linestyle='-', color='r',markerfacecolor='b')

# Adding labels and title
plt.xlabel('x')
plt.ylabel('P(x)')
plt.title('P(x) vs x')


# Display the plot
plt.show()


# Question (b)

# Approximated roots list is below:
roots_lst = [0.035,0.165,0.385,0.615,0.835,0.965]
exct_roots = []
for i in roots_lst:
    x = i
    while abs(P(x)) > 10**-10:
        x = x - P(x)/deriv_P(x)
        
    exct_roots.append(x)


print("The roots are")
for i in exct_roots:
    print("{:4.12f}".format(i))

print("complete")

# To see the values of roots after running the code, close the graph plotted.
             