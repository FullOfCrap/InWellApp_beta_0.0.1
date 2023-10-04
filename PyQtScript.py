import math

num1 = 20
num2 = 1.35
num3 = 0.350
num4 = 20
num6 = 11
num7 = 18
num8 = 100
num5 = 20
num9 = float(num3) / 2
b = 0.86

K = (float(num1) * (math.log10(float(num8) / float(num9)))) / (
        (1.366 * ((2 * float(num6)) - float(num2))) * float(num2)) * (1 / b)

print(b)
print(K)
K = K / 3600
R = 3000 * float(num2) * math.sqrt(K)
print(K)
print(R)
K = (float(num1) * (math.log10(float(num8) / float(num9)))) / (
        (1.366 * ((2 * float(num6)) - float(num2))) * float(num2)) * (1 / b)
K = K / 3600
R = 3000 * float(num2) * math.sqrt(K)
print(K)
print(R)
