import math
import csv
import random
import time
import numpy as np
import numpy as p
import pandas as pd
import matplotlib.pyplot as plt
from numpy import linalg as LA
from sympy import simplify, Symbol, exp


def sum_x(points_f):
    sum_4 = sum_3 = sum_2 = sum_1 = 0
    for _ in range(len(points_f)):
        sum_4 += points_f[_][0] ** 4
        sum_3 += points_f[_][0] ** 3
        sum_2 += points_f[_][0] ** 2
        sum_1 += points_f[_][0]
    return sum_4, sum_3, sum_2, sum_1


def sum_y(points_f):
    sum_y2 = sum_y1 = y_sum = 0
    for _ in range(len(points_f)):
        sum_y2 += points_f[_][0] ** 2 * points_f[_][1]
        sum_y1 += points_f[_][0] * points_f[_][1]
        y_sum += points_f[_][1]
    return sum_y2, sum_y1, y_sum


def gauss1(matrix):
    """
    Функция решает матрицу методом Гаусса

    :param matrix: list -- матрица
    :return matrix: решенная матрица
    """
    n_f = len(matrix)
    for k_f in range(n_f):
        for _ in range(k_f, n_f):
            if _ == k_f:
                matrix[k_f][len(matrix[_]) - 1] = matrix[k_f][
                                                      len(matrix[_]) - 1] / \
                                                  matrix[k_f][k_f]
            else:
                matrix[_][len(matrix[_]) - 1] = matrix[_][len(matrix[_]) - 1]
                - ((matrix[k_f][len(matrix[_]) - 1]) * (matrix[_][k_f]))

            for index_f in range(len(matrix[_]) - 2, -1, -1):
                if _ == k_f:
                    matrix[_][index_f] = matrix[_][index_f] / matrix[_][_]
                else:
                    matrix[_][index_f] = \
                        matrix[_][index_f] - \
                        matrix[k_f][index_f] * matrix[_][k_f]
    return matrix


def kvadr(points, n):
    """

    :param points: list -- значения x и y
    :param n: int -- количество точек
    :return x_f: параметры квадратичной функции
    """
    matrix = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    x_f = list(sum_x(points))
    y = sum_y(points)
    x_f.append(n)
    for _ in range(3):
        q = _
        for index_f in range(3):
            matrix[_][index_f] = x_f[q]
            q += 1
        matrix[_][index_f + 1] = y[_]
    matrix2 = gauss1(matrix)
    matrix_a = []
    matrix_b = []
    for _ in range(3):
        c = []
        b = []
        for index_f in range(3):
            c.append(matrix2[_][index_f])
            b.append(matrix[_][index_f])
        matrix_a.append(c)
        matrix_b.append(matrix2[_][index_f + 1])
    matrix_a_obr = LA.inv(matrix_a)
    x_f = np.dot(matrix_a_obr, matrix_b)
    return x_f


def type_number(ch):
    if '.' in ch:
        ch = float(ch)
    else:
        ch = int(ch)
    return ch


def ravn(sp):
    """
    Функция вычисляет равностоящие узлы интерполяции или нет
    :param sp: list -- знчения x и y
    :return p_f: число, которое означает равностоящие узлы интерполяции или нет
    """
    f = sp[0][0] - sp[1][0]
    p_f = 1
    for _ in range(1, len(sp) - 1):
        if (sp[_][0] - sp[_ + 1][0]) != f:
            p_f = 0
            break
    return p_f


def fun(points):
    """
    Функцией нормального распределения
    :param points: list -- значения x и y
    :return: параметры функции
    """
    x_ = sum([points[_][0] for _ in range(len(points))]) / len(points)
    s_ = 0
    b = 0
    for _ in range(len(points)):
        s_ += (points[_][0] - x_) ** 2
        b += points[_][0] / len(points)
    c = math.sqrt(s_ / len(points))
    a = 1 / (c * math.sqrt(2 * math.pi))
    return a, b, c


def calc_a(k):
    return delta[k][0] / (math.factorial(k) * (h ** k))


def Lagrang(points, r, x):
    """
    Нахождение значения функции в точке x
    интерполяционным многочленом Лагранжа

    :param points: list -- значения x и y
    :param r: int -- число, которое означает равностоящие узлы интерполяции или нет
    :param x: float -- точка, в которой нужно найти значение функции
    :return L: значение функции в точке x
    """
    L = 0
    n = len(points)
    if r == 0:
        for _ in range(n):
            L += Lagrang1(x, points[_][0], points[_][1], _, points)

    else:
        q = (x - points[0][0]) / (points[0][0] - points[1][0])
        for _ in range(n):
            L += Lagrang2(q, points[_][1], _, n)
        for _ in range(n + 1):
            L *= (q - _)
        L = L / math.factorial(n)
    return L


def Lagrang1(x, xi, yi, ind, t):
    """
    Вычисление множителей многочлена Лагранжа
    если узлы интерполяции неравноотстоящие

    :param x: float -- точка, в которой нужно найти значение функции
    :param xi: float -- значение x
    :param yi: float -- значение y
    :param ind: int -- индекс точки в списке
    :param t: list -- значения x и y
    :return Ln: множитель многочлена Лагранжа
    """
    Ln = yi
    for e in range(len(t)):
        if e != ind:
            if (xi - t[e][0]) == 0:
                Ln *= 0
            else:
                Ln *= ((x - t[e][0]) / (xi - t[e][0]))
    return Ln


def Lagrang2(q, yi, ind, n):
    """
    Вычисление множителей многочлена Лагранжа
    если узлы интерполяции равноотстоящие

    :param q: float -- множитель многочлена Лагранжа, вычесленный по формуле
    :param yi: значение y
    :param ind: int -- индекс точки в списке
    :param n: int -- количество точек
    :return Ln: множитель многочлена Лагранжа
    """
    Ln = (-1) ** (n - ind) * yi * (math.factorial(n) / (
            math.factorial(ind) * math.factorial(n - ind) * (q - ind)))
    return Ln


def approx1(m1):
    """
    Выводит значения yi аппроксимированной функции

    :param m1: list -- значения x и y
    :return f1: значения yi аппроксимированной функции
    """

    s1 = 0
    s2 = 0
    s3 = 0
    s4 = 0
    s5 = 0
    for k in range(len(m1)):
        s1 += m1[k, 0]
        s2 += m1[k, 1]
        s3 += m1[k, 0] ** 2
        s4 += m1[k, 0] * m1[k, 1]
        s5 += 1
    det = LA.det(np.array([[s3, s1], [s1, s5]]))
    det1 = LA.det(np.array([[s4, s1], [s2, s5]]))
    det2 = LA.det(np.array([[s3, s4], [s1, s2]]))
    a = det1 / det
    b = det2 / det
    fi = []
    f1 = []
    for _ in range(len(m1)):
        fi.append([m1[_, 0], m1[_, 1], a * m1[_, 0] + b])
        f1.append(a * m1[_, 0] + b)
    sdis = 0
    for _ in range(len(m1)):
        sdis += (m1[_, 1] - f1[_]) ** 2
    return f1


def xi(m1):
    """
    Значения xi
    :param m1: list -- значения x и y
    :return mas: значение xi
    """
    mas = []
    for _ in range(len(m1)):
        mas.append(m1[_, 0])
    return mas


def yi(m1):
    """
    Значения yi
    :param m1: list -- значения x и y
    :return mas1: значение yi
    """
    mas1 = []
    for i in range(len(m1)):
        mas1.append(m1[i, 1])
    return mas1


def inter(m1):
    """
    Интерполяция Методом Ньютона
    :param m1: list -- значения x и y
    :return f2: значения функции
    """
    h = m1[1, 0] - m1[0, 0]
    n = len(m1)
    delta = []
    po = []
    for p in range(len(m1)):
        po.append(m1[p, 1])
    delta.append(po)
    delta = list(delta)
    for _ in range(1, n + 1):
        current = []
        for j in range(n - _):
            current.append(delta[_ - 1][j + 1] - delta[_ - 1][j])
        delta.append(current)
    delta.append([0])

    def calc_a(k):
        return delta[k][0] / (math.factorial(k) * (h ** k))

    polynom = np.poly1d([0])
    for k in range(n):
        a = calc_a(k)
        cur_polynom = np.poly1d([1])
        for j in range(k):
            cur_polynom *= np.poly1d([1, -m1[j, 0]])
        cur_polynom = a * cur_polynom
        polynom += cur_polynom

    fi1 = []
    f2 = []
    for _ in range(len(m1)):
        fi1.append([m1[_, 0], m1[_, 1], polynom(m1[_, 0])])
        f2.append(polynom(m1[_, 0]))
    return f2


def inter1(m1):
    """
    Интерполяция Методом Ньютона назад
    :param m1: list -- значения x и y
    :return f2: значения функции
    """
    h = m1[1, 0] - m1[0, 0]
    n = len(m1)
    delta = []
    po = []
    for p in range(len(m1)):
        po.append(m1[p, 1])
    delta.append(po)

    for _ in range(1, n + 1):
        current = []
        for index in range(n - _):
            current.append(delta[_ - 1][index + 1] - delta[_ - 1][index])
        delta.append(current)
    delta.append([0])

    def calc_a(k):
        return delta[n - 1 - k][-1] / (
                math.factorial(n - 1 - k) * (h) ** (n - 1 - k))

    polynom = np.poly1d([0])
    for k in range(n - 1, -1, -1):
        a = calc_a(k)
        cur_polynom = np.poly1d([1])
        for index in range(n - 1, k - 1, -1):
            cur_polynom *= np.poly1d([1, -m1[index, 0]])
        cur_polynom = a * cur_polynom
        polynom += cur_polynom

    fi1 = []
    f2 = []
    for _ in range(len(m1)):
        fi1.append([m1[_, 0], m1[_, 1], polynom(m1[_, 0])])
        f2.append(polynom(m1[_, 0]))
    return f2


def interp(x, xp, yp):
    """
    Интерполяция функцией библиотеки numpy
    :param x: float -- значение x
    :param xp: list -- значения x
    :param yp: list -- значения y
    :return: значения функции
    """
    return np.interp(x, xp, yp)


def numpy_approximation(x_values, y_values):
    """
    Аппроксимация функцией библиотеки numpy
    :param x_values: list -- значения x
    :param y_values: list -- значения y
    :return: значения аппроксимирующей функции
    """

    def genetate_poly(coef):
        coef.reverse()

        def poly(x):
            return sum(coef[i] * x ** i for i in range(len(coef)))

        return poly

    f = genetate_poly(list(
        np.polynomial.polynomial.Polynomial.fit(x_values, y_values, 10).coef))
    f_values = [f(x) for x in x_values]
    return x_values, y_values, f_values


def read_excel(name):
    """
    Считывание CSV айла
    :param name: str -- название файла
    :return points: значения x и y из файла
    """
    curs = pd.read_excel(name)
    y = curs.curs.values.tolist()
    x = [i for i in range(len(y))]
    points = []
    for i in range(len(x)):
        points.append([x[i], y[i]])
    return points


n = int(input('Введите количество точек '))
print('Введите диапазон, из которого будут генерироваться числа')
a = float(input('a = '))
b = float(input('b = '))
print("Какие числа хотите генерировать?")
print("1-Целые")
print("2-Дробные")
t = int(input('Введите число '))
points = []
x = float(input('Введите х '))

with open("matrix.csv", mode="w", encoding='utf-8') as file:
    file_writer = csv.writer(file, delimiter=",", lineterminator="\r")
    for i in range(n):
        st = []
        for j in range(2):
            if t == 1:
                st.append(int(random.uniform(a, b)))
            else:
                st.append(float(random.uniform(a, b)))
        file_writer.writerow(st)

with open("matrix.csv") as file:
    reader = csv.reader(file, delimiter=',', quotechar=',')
    for row in reader:
        points.append(row)

for _ in range(n):
    for index in range(2):
        points[_][index] = type_number(points[_][index])

points.sort()
r = ravn(points)
point = []
for _ in range(n):
    st = []
    for j in range(2):
        st.append(Symbol(str(points[_][j])))
    point.append(st)

L = Lagrang(points, r, x)
Lstr = Lagrang(point, r, (Symbol(str(x))))

f = []
for _ in range(n):
    f.append(Lagrang(points, r, points[_][0]))
    points[_].append(f[_])
print('Массив  точек')
print(points)
print(L)
print('Интерполяционный многочлен Лагранжа')
Lstr

pl = plt.figure(figsize=(14, 11))
plt.plot([points[i][0] for i in range(n)], [points[i][1] for i in range(n)],
         'ro')
plt.plot([points[i][0] for i in range(n)], f, 'purple')
plt.legend(['узлы', 'f(x)'], fontsize=13, loc=2)
plt.xlabel("Ось X", fontsize=14, color='black')
plt.ylabel("Ось Y", fontsize=14, color='black')
plt.show()

c3 = int(input("Введите значение, от которого будет изменятся x"))
c4 = int(input("Введите значение, до которого будет изменятся x"))
h1 = int(input("Введите шаг по x"))

with open("mrx.csv", mode="w", encoding='utf-8') as file:
    file_writer = csv.writer(file, delimiter=",", lineterminator="\r")

    a1 = []
    for i in range(c3, c4, h1):
        a1.append(i)
    file_writer.writerow(a1)
    a1 = []
    for j in range(c3, c4, h1):
        a1.append(int(random.uniform(0, 100)))
    file_writer.writerow(a1)

m1 = []

with open("mrx.csv") as file:
    reader = csv.reader(file, delimiter=',')
    for row in reader:
        m1.append(row)

m1 = np.array(m1, np.int32)
print(m1)

h = m1[0, 1] - m1[0, 0]
n = len(m1[0])
delta = []
delta.append(list(m1[1]))
for i in range(1, n + 1):
    current = []
    for j in range(n - i):
        current.append(delta[i - 1][j + 1] - delta[i - 1][j])
    delta.append(current)
delta.append([0])

polynom = np.poly1d([0])
for k in range(n):
    a = calc_a(k)
    cur_polynom = np.poly1d([1])
    for j in range(k):
        cur_polynom *= np.poly1d([1, -m1[0, j]])
    cur_polynom = a * cur_polynom
    polynom += cur_polynom

print("Интерполяционный многочлен Ньютона")
print(polynom)
fi1 = []
f2 = []
for i in range(len(m1[0])):
    fi1.append([m1[0, i], m1[1, i], polynom(m1[0, i])])
    f2.append(polynom(m1[0, i]))
print("Массив точек в формате [xi,yi,fi]")
print(fi1)


plt.plot(m1[0], m1[1], 'r', label='f*(x) первон.')
plt.plot(m1[0], f2, 'm', label='f*(x) интер.')
plt.show()

c1 = int(input("Введите значение, от которого будет изменятся x"))
c2 = int(input("Введите значение, до которого будет изменятся x"))
h = int(input("Введите шаг по x"))

with open("mrx.csv", mode="w", encoding='utf-8') as file:
    file_writer = csv.writer(file, delimiter=",", lineterminator="\r")

    a1 = []
    for i in range(c1, c2, h):
        a1.append(i)
    file_writer.writerow(a1)
    a1 = []
    for j in range(c1, c2, h):
        a1.append(int(random.uniform(0, 100)))
    file_writer.writerow(a1)

m1 = []

# открываем файл и записываем в матрицу
with open("mrx.csv") as file:
    reader = csv.reader(file, delimiter=',')
    for row in reader:
        m1.append(row)
m1 = np.array(m1, np.int32)
print(m1)

s1 = 0
s2 = 0
s3 = 0
s4 = 0
s5 = 0
for k in range(len(m1[0])):
    s1 += m1[0, k]
    s2 += m1[1, k]
    s3 += m1[0, k] ** 2
    s4 += m1[0, k] * m1[1, k]
    s5 += 1

det = LA.det(np.array([[s3, s1], [s1, s5]]))
det1 = LA.det(np.array([[s4, s1], [s2, s5]]))
det2 = LA.det(np.array([[s3, s4], [s1, s2]]))
a = det1 / det
b = det2 / det

print("Вид аппроксимирующей функции")
if b > 0:
    print('y=', a, '*x+', b, )
else:
    print('y=', a, '*x', b, )

fi = []
f1 = []
for i in range(len(m1[0])):
    fi.append([m1[0, i], m1[1, i], a * m1[0, i] + b])
    f1.append(a * m1[0, i] + b)
print("Массив точек в формате [xi,yi,fi]")
print(fi)

sdis = 0
for i in range(len(m1[0])):
    sdis += (m1[1, i] - f1[i]) ** 2
print("Дисперсия")
print(sdis)

plt.plot(m1[0], m1[1], 'r', label='f*(x) первон.')
plt.plot(m1[0], f1, 'm', label='f*(x) аппрокс.')
plt.show()

n = int(input('Введите количество точек '))
print('Введите диапазон, из которого будут генерироваться числа')
a = float(input('a = '))
b = float(input('b = '))
print("Какие числа хотите генерировать?")
print("1-Целые")
print("2-Дробные")
t = int(input('Введите число '))
points = []

with open("matrix.csv", mode="w", encoding='utf-8') as file:
    file_writer = csv.writer(file, delimiter=",", lineterminator="\r")
    for _ in range(n):
        st = []
        for j in range(2):
            if t == 1:
                st.append(int(random.uniform(a, b)))
            else:
                st.append(float(random.uniform(a, b)))
        file_writer.writerow(st)

with open("matrix.csv") as file:
    reader = csv.reader(file, delimiter=',', quotechar=',')
    for row in reader:
        points.append(row)

for _ in range(n):
    for j in range(2):
        points[_][j] = type_number(points[_][j])
points.sort()
X = kvadr(points, n)
Xs = []
f = []
x_1 = []
d = 0
for i in range(n):
    f.append(X[0] * points[i][0] ** 2 + X[1] * points[i][0] + X[2])
    points[i].append(f[i])
    d += (points[i][1] - f[i]) ** 2
for i in range(len(X)):
    Xs.append(Symbol(str(X[i])))
print('Массив  точек')
print(points)
print('Величина дисперсии')
print(d)
print('Вид аппроксимирующей функции')
simplify(Xs[0] * Symbol('x') ** 2 + Xs[1] * Symbol('x') + Xs[2])

plt.plot([points[i][0] for i in range(n)], [points[i][1] for i in range(n)],
         'blue')
plt.plot([points[i][0] for i in range(n)], f, 'purple')
plt.legend(['f*(x) первон.', 'f*(x) аппрокс.'], fontsize=13, loc=2)
plt.show()

n = int(input('Введите количество точек '))
print('Введите диапазон, из которого будут генерироваться числа')
a = float(input('a = '))
b = float(input('b = '))
print("Какие числа хотите генерировать?")
print("1-Целые")
print("2-Дробные")
t = int(input('Введите число '))
points = []
with open("matrix.csv", mode="w", encoding='utf-8') as file:
    file_writer = csv.writer(file, delimiter=",", lineterminator="\r")
    for i in range(n):
        st = []
        for j in range(2):
            if t == 1:
                st.append(int(random.uniform(a, b)))
            else:
                st.append(float(random.uniform(a, b)))
        file_writer.writerow(st)

with open("matrix.csv") as file:
    reader = csv.reader(file, delimiter=',', quotechar=',')
    for row in reader:
        points.append(row)

for _ in range(n):
    for j in range(2):
        points[_][j] = type_number(points[_][j])
points.sort()
X = fun(points)
f = []
d = 0
for _ in range(n):
    f.append(X[0] * exp(-(points[_][0] - X[1]) ** 2 / (X[2]) ** 2))
    points[_].append(f[_])
    d += (points[_][1] - f[_]) ** 2
print('Аппроксимация функцией нормального распределения:')
print('Массив  точек ')
print(points)
print('Величина дисперсии ')
print(d)
print('Вид аппроксимирующей функции  ')
simplify(X[0] * exp(-(Symbol('x') - X[1]) ** 2 / (X[2]) ** 2))

plt.plot([points[_][0] for _ in range(n)], [points[_][1] for _ in range(n)],
         'blue')
plt.plot([points[_][0] for _ in range(n)], f, 'purple')
plt.legend(['f*(x) первон.', 'f*(x) аппрокс.'], fontsize=13, loc=2)
plt.show()

m1 = [[1, 10], [2, 55], [3, 71]]
m1 = np.array(m1)

n = int(input('Введите количество точек '))
print('Введите диапазон, из которого будут генерироваться числа')
a = float(input('a = '))
b = float(input('b = '))
print("Какие числа хотите генерировать?")
print("1-Целые")
print("2-Дробные")
t = int(input('Введите число '))
points = []
with open("matrix.csv", mode="w", encoding='utf-8') as file:
    file_writer = csv.writer(file, delimiter=",", lineterminator="\r")
    for i in range(n):
        st = []
        for j in range(2):
            if t == 1:
                st.append(int(random.uniform(a, b)))
            else:
                st.append(float(random.uniform(a, b)))
        file_writer.writerow(st)

with open("matrix.csv") as file:
    reader = csv.reader(file, delimiter=',', quotechar=',')
    for row in reader:
        points.append(row)

for i in range(n):
    for j in range(2):
        points[i][j] = type_number(points[i][j])
points.sort()


plt.plot(xi(np.array(points)), yi(np.array(points)), 'pink',
         label='f*(x) первон.')
plt.plot(xi(np.array(points)), approx1(np.array(points)), 'purple',
         label='f*(x) аппрокс.')
plt.title('Аппроксимация линейной функцией', fontsize=10,
          fontname='Times New Roman')
plt.legend(['f*(x) первон.', 'f*(x) аппрокс.'], fontsize=9, loc=2)
plt.xlabel("Ось X", fontsize=14, color='black')
plt.ylabel("Ось Y", fontsize=14, color='black')
plt.show()

X = kvadr(points, n)
f = []
x_1 = []
d = 0

for i in range(n):
    f.append(X[0] * points[i][0] ** 2 + X[1] * points[i][0] + X[2])
    points[i].append(f[i])
    d += (points[i][1] - f[i]) ** 2

print('Величина дисперсии')
print(d)

plt.plot([points[i][0] for i in range(n)], [points[i][1] for i in range(n)],
         'pink')
plt.plot([points[i][0] for i in range(n)], f, 'purple')
plt.title('Аппроксимация квадратичной функцией', fontsize=10,
          fontname='Times New Roman')
plt.legend(['f*(x) первон.', 'f*(x) аппрокс.'], fontsize=9, loc=2)
plt.xlabel("Ось X", fontsize=14, color='black')
plt.ylabel("Ось Y", fontsize=14, color='black')
plt.show()

X = fun(points)
f = []
d = 0
for i in range(n):
    f.append(X[0] * math.exp(-(points[i][0] - X[1]) ** 2 / (X[2]) ** 2))
    points[i].append(f[i])
    d += (points[i][1] - f[i]) ** 2

plt.plot([points[i][0] for i in range(n)], [points[i][1] for i in range(n)],
         'pink')
plt.plot([points[i][0] for i in range(n)], f, 'purple')
plt.title('Аппроксимация функцией нормального распределения', fontsize=10,
          fontname='Times New Roman')
plt.legend(['f*(x) первон.', 'f*(x) аппрокс.'], fontsize=9, loc=2)
plt.xlabel("Ось X", fontsize=14, color='black')
plt.ylabel("Ось Y", fontsize=14, color='black')
plt.show()
print('Величина дисперсии')
print(d)

A = numpy_approximation([points[i][0] for i in range(n)],
                        [points[i][1] for i in range(n)])
plt.plot(A[0], A[1], 'pink')
plt.plot(A[0], A[2], 'purple')
plt.title('Аппроксимация функцией библиотеки numpy', fontsize=10,
          fontname='Times New Roman')
plt.legend(['f*(x) первон.', 'f*(x) аппрокс.'], fontsize=9, loc=2)
plt.xlabel("Ось X", color='black')
plt.ylabel("Ось Y", color='black')
plt.show()

f = []
for i in range(n):
    f.append(Lagrang(points, r, points[i][0]))
    points[i].append(f[i])

plt.plot([points[i][0] for i in range(n)], [points[i][1] for i in range(n)],
         'ro')
plt.plot([points[i][0] for i in range(n)], f, 'blue')
plt.title('Интерполяция методом Лагранжа', fontsize=10,
          fontname='Times New Roman')
plt.legend(['f*(x) первон.', 'f*(x) интер.'], fontsize=9, loc=2)
plt.xlabel("Ось X", fontsize=14, color='black')
plt.ylabel("Ось Y", fontsize=14, color='black')
plt.show()

plt.plot(xi(np.array(points)), yi(np.array(points)), 'blue',
         label='f*(x) первон.')
plt.plot(xi(np.array(points)), inter(np.array(points)), 'orange',
         label='f*(x) аппрокс.')
plt.title('Интерполяция методом Ньютона вперед', fontsize=10,
          fontname='Times New Roman')
plt.legend(['f*(x) первон.', 'f*(x) интер.'], fontsize=9, loc=2)
plt.xlabel("Ось X", fontsize=14, color='black')
plt.ylabel("Ось Y", fontsize=14, color='black')
plt.show()


plt.plot(xi(np.array(points)), yi(np.array(points)), 'blue',
         label='f*(x) первон.')
plt.plot(xi(np.array(points)), inter1(np.array(points)), 'orange',
         label='f*(x) аппрокс.')
plt.title('Интерполяция методом Ньютона назад', fontsize=10,
          fontname='Times New Roman')
plt.legend(['f*(x) первон.', 'f*(x) интер.'], fontsize=9, loc=2)
plt.xlabel("Ось X", fontsize=14, color='black')
plt.ylabel("Ось Y", fontsize=14, color='black')
plt.show()

xp = [points[i][0] for i in range(n)]
yp = [points[i][1] for i in range(n)]
f = []

for i in range(n):
    f.append([points[i][0], points[i][1], interp(points[i][0], xp, yp)])

plt.plot([points[i][0] for i in range(n)], [points[i][1] for i in range(n)],
         'ro')
plt.plot([points[i][0] for i in range(n)], [f[i][2] for i in range(n)], 'blue')
plt.title('Интерполяция функцией библиотеки numpy', fontsize=10,
          fontname='Times New Roman')
plt.legend(['f*(x) первон.', 'f*(x) интер.'], fontsize=9, loc=2)
plt.show()


def read_excel(name):
    curs = pd.read_excel(name)
    y = curs.curs.values.tolist()
    x_f = [i for _ in range(len(y))]
    points_f = []
    for _ in range(len(x_f)):
        points_f.append([x_f[_], y[_]])
    return points_f


points = read_excel('/Users/tkach/Downloads/валюта.xlsx')
n = (len(points))

X = fun(points)
f = []
d = 0
for i in range(n):
    f.append(X[0] * exp(-(points[i][0] - X[1]) ** 2 / (X[2]) ** 2))
    points[i].append(f[i])
    d += (points[i][1] - f[i]) ** 2

plt.plot([points[i][0] for i in range(n)], [points[i][1] for i in range(n)],
         'blue')
plt.plot([points[i][0] for i in range(n)], f, 'purple')
plt.title('Аппроксимация функцией нормального распределения', fontsize=10,
          fontname='Times New Roman')
plt.legend(['курс доллара С 01.12.2013 по 01.01.2015', 'f*(x) аппрокс.'],
           fontsize=9, loc=2)
plt.show()
