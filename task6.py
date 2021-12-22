import numpy as np
import sympy
from sympy import *
import sympy.abc
from sympy.abc import *
from sympy import symbols
import re
import matplotlib.pyplot as plt
import pandas as pd
from numpy import linalg as LA
import math


def sum_x(tochki):
    S4 = S3 = S2 = S1 = 0
    for i in range(len(tochki)):
        S4 += tochki[i][0] ** 4
        S3 += tochki[i][0] ** 3
        S2 += tochki[i][0] ** 2
        S1 += tochki[i][0]
    return S4, S3, S2, S1


def sum_y(tochki):
    Sy2 = Sy1 = Sy = 0
    for i in range(len(tochki)):
        Sy2 += tochki[i][0] ** 2 * tochki[i][1]
        Sy1 += tochki[i][0] * tochki[i][1]
        Sy += tochki[i][1]
    return Sy2, Sy1, Sy


def gauss1(matrix):
    n = len(matrix)
    for k in range(n):
        for i in range(k, n):
            if i == k:
                matrix[k][len(matrix[i]) - 1] =\
                    matrix[k][len(matrix[i]) - 1] / matrix[k][k]
            else:
                matrix[i][len(matrix[i]) - 1] = \
                    matrix[i][len(matrix[i]) - 1] - \
                    ((matrix[k][len(matrix[i]) - 1]) * (matrix[i][k]))

            for j in range(len(matrix[i]) - 2, -1, -1):
                if i == k:
                    matrix[i][j] = matrix[i][j] / matrix[i][i]
                else:
                    matrix[i][j] = matrix[i][j] - matrix[k][j] * matrix[i][k]
    return matrix


def kvadr_diff(tochki, n):
    """
    Функция аппроксимации
    """
    matrix = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    x = list(sum_x(tochki))
    y = sum_y(tochki)
    x.append(n)
    for i in range(3):
        count = i
        for j in range(3):
            matrix[i][j] = x[count]
            count += 1
        matrix[i][j + 1] = y[i]
    matrix2 = gauss1(matrix)
    matrix_A = []
    matrix_B = []
    for i in range(3):
        c = []
        b = []
        for j in range(3):
            c.append(matrix2[i][j])
            b.append(matrix[i][j])
        matrix_A.append(c)
        matrix_B.append(matrix2[i][j + 1])
    matrix_A_obr = LA.inv(matrix_A)
    coord = np.dot(matrix_A_obr, matrix_B)
    return (coord)


def inter_diff(m1):
    """
    Функция интерполяции методом Ньютона
    """
    h = m1[1, 0] - m1[0, 0]
    n = len(m1)
    delta = []
    po = []
    for p in range(len(m1)):
        po.append(m1[p, 1])
    delta.append(po)
    delta = list(delta)
    for i in range(1, n + 1):
        current = []
        for j in range(n - i):
            current.append(delta[i - 1][j + 1] - delta[i - 1][j])
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
    for i in range(len(m1)):
        fi1.append([m1[i, 0], m1[i, 1], polynom(m1[i, 0])])
        f2.append(polynom(m1[i, 0]))
    return (f2)


def approximation_diff(b, variable_new, name):
    """
    Функция для построения графиков аппроксимации
    :param b: list -- список с решением уравения вида: [x,y]
    :param variable_new: list -- список переменных уравнения
    :param name: str -- название метода

    :return x_list: list -- список интерполированных точек
    """
    print("Аппроксимация ", name)
    for j in range(len(b[0]) - 2):
        tochki = []
        for i in range(len(b)):
            qp = []
            qp.append(b[i][1])
            qp.append(b[i][j + 2])
            tochki.append(qp)
        n = len(tochki)
        x_list = kvadr_diff(tochki, n)
        fx = []
        x_1 = []
        for ii in range(n):
            fx.append(
                x_list[0] * tochki[ii][0] ** 2 + x_list[1] * tochki[ii][0] +
                x_list[2])
            tochki[ii].append(fx[ii])
        plt.title(variable_new[j], fontsize=18, fontname='Times New Roman')
        plt.plot([tochki[i][0] for i in range(n)],
                 [tochki[i][1] for i in range(n)],
                 'pink')
        plt.plot([tochki[i][0] for i in range(n)], fx, 'blue')
        plt.legend(['f*(x) первон.', 'f*(x) аппрокс.'], fontsize=13, loc=2)
        plt.show()
    return x_list


def interpol_diff(b, variable_new, name):
    """
    Функция для построения графиков интерполяции
    :param b: list -- список с решением уравения вида: [x,y]
    :param variable_new: list -- список переменных уравнения
    :param name: str -- название метода

    :return x_list: list -- список интерполированных точек
    """
    print("Интерполяция методом Ньютона ", name)
    for j in range(len(b[0]) - 2):
        tochki = []
        for i in range(len(b)):
            qp = []
            qp.append(b[i][1])
            qp.append(b[i][j + 2])
            tochki.append(qp)
        x_list = inter_diff(np.array(tochki))
        n = len(tochki)
        plt.title(variable_new[j], fontsize=18, fontname='Times New Roman')
        plt.plot([tochki[i][0] for i in range(n)],
                 [tochki[i][1] for i in range(n)],
                 'pink')
        plt.plot([tochki[i][0] for i in range(n)], x_list, 'blue')
        plt.legend(['f*(x) первон.', 'f*(x) аппрокс.'], fontsize=13, loc=2)
        plt.show()
    return x_list


def Euler(variable_new, x0, start, n, ab, f):
    """
    Функция для решения дифф. ур. методом Эйлера

    :param variable_new: list -- список переменных уравения
    :param x0: float -- начальное условие x
    :param start: list -- начальные условия переменных в точке x0
    :param n: int -- количество уравнений
    :param ab: list -- промежуток x
    :param f: function _lambdifygenerated --
    список функций, выполняющих подстановку значений в уравнения

    :return solution: list -- список - таблица с решением ур вида:
    x0, y, y', ..., z, z' ...
    """
    h = (ab[1] - ab[0]) / n
    solution = [[0 for i in range(len(variable_new))] for j in range(n + 1)]
    for i in range(len(variable_new)):
        solution[0][i] = start[i]
    xl = [0] * (n + 1)
    xl[0] = x0
    for i in range(1, n + 1):
        for j in range(len(variable_new)):
            step = h * f[j](xl[i - 1], solution[i - 1])
            solution[i][j] = solution[i - 1][j] + step
        xl[i] = xl[i - 1] + h
    for i in range(n + 1):
        solution[i].insert(0, xl[i])
        solution[i].insert(0, i)
    return solution


def Euler_Kauchi(variable_new, x0, start, n, ab, f):
    """
   Функция для решения дифф. ур. методом Эйлера

    :param variable_new: list -- список переменных уравения
    :param x0: float -- начальное условие x
    :param start: list -- начальные условия переменных в точке x0
    :param n: int -- количество уравнений
    :param ab: list -- промежуток x
    :param f: function _lambdifygenerated --
    список функций, выполняющих подстановку значений в уравнения

    :return solution: list -- список - таблица с решением ур вида:
    x0, y, y', ..., z, z' ...
    """
    h = (ab[1] - ab[0]) / n
    solution = [[0 for i in range(len(variable_new))] for j in range(n + 1)]
    for i in range(len(variable_new)):
        solution[0][i] = start[i]
    xl = [0] * (n + 1)
    xl[0] = x0
    for i in range(1, n + 1):
        f0 = []
        for g in range(0, len(variable_new)):
            f0.append(solution[i - 1][g] + h * f[g](xl[i - 1],
                                                    solution[i - 1]))
        xl[i] = xl[i - 1] + h
        for j in range(len(variable_new)):
            step = 0.5 * h *\
                   (f[j](xl[i - 1], solution[i - 1]) + f[j](xl[i], f0))
            solution[i][j] = solution[i - 1][j] + step
    for i in range(n + 1):
        solution[i].insert(0, xl[i])
        solution[i].insert(0, i)
    return solution


def koef(f, xl, solution, h, variable_new):
    """
    k1, k2, k3, k4 для Метода Рунгу-Кутты
    """
    k = [[0 for i in range(len(variable_new))] for j in range(4)]
    # k1
    for j in range(0, len(variable_new)):
        k[0][j] = h * f[j](xl, solution)
    # k2
    for j in range(0, len(variable_new)):
        k[1][j] = h * f[j](xl + h / 2,
                           list(np.array(solution) + np.array(k[0]) / 2))
    # k3
    for j in range(0, len(variable_new)):
        k[2][j] = h * \
                  f[j](xl + h / 2,
                       list(np.array(solution) + np.array(k[1]) / 2))
    # k4
    for j in range(0, len(variable_new)):
        k[3][j] = h * f[j](xl + h,
                           list(np.array(solution) + np.array(k[2])))
    return k


def Runge_Kutta(variable_new, x0, start, n, ab, f):
    """"
    Функция для решения дифф. ур. методом Эйлера

    :param variable_new: list -- список переменных уравения
    :param x0: float -- начальное условие x
    :param start: list -- начальные условия переменных в точке x0
    :param n: int -- количество уравнений
    :param ab: list -- промежуток x
    :param f: function _lambdifygenerated --
    список функций, выполняющих подстановку значений в уравнения

    :return solution: list -- список - таблица с решением ур вида:
    [x0, y, y', ..., z, z' ...]
    """
    h = (ab[1] - ab[0]) / n
    solution = [[0 for i in range(len(variable_new))] for j in range(n + 1)]
    for i in range(len(variable_new)):
        solution[0][i] = start[i]
    xl = [0] * (n + 1)
    xl[0] = x0
    for i in range(1, n + 1):
        kf = koef(f, xl[i - 1], solution[i - 1], h, variable_new)
        for j in range(0, len(variable_new)):
            solution[i][j] =\
                solution[i - 1][j] + (1 / 6) * \
                (kf[0][j] + 2 * kf[1][j] + 2 * kf[2][j] + kf[3][j])
        xl[i] = xl[i - 1] + h
    for i in range(n + 1):
        solution[i].insert(0, xl[i])
        solution[i].insert(0, i)
    return solution


def equation_2(euler, euler_cauchy, runge_kutta):
    """
    Функция для построения графиков для 2 ур

    :param euler: list -- список с решением уравнения методом Эйлера
    вида: [x,y,z]
    :param euler_cauchy: list -- список с решением уравнения
    методом Эйлера-Коши вида: [x,y,z]
    :param runge_kutta: list -- список с решением уравнения
    методом Рунге-Кутты вида: [x,y,z]

    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    x = [euler[i][1] for i in range(len(euler))]
    y = [euler[i][2] for i in range(len(euler))]
    z = [euler[i][3] for i in range(len(euler))]
    ax.plot(x, y, z, label='parametric curve', color='purple', linewidth=3)
    x = [euler_cauchy[i][1] for i in range(len(euler_cauchy))]
    y = [euler_cauchy[i][2] for i in range(len(euler_cauchy))]
    z = [euler_cauchy[i][3] for i in range(len(euler_cauchy))]
    ax.plot(x, y, z,
            label='parametric curve',
            marker='X', markersize=5,
            color='teal',
            markerfacecolor='greenyellow',
            markeredgecolor='teal')
    x = [runge_kutta[i][1] for i in range(len(runge_kutta))]
    y = [runge_kutta[i][2] for i in range(len(runge_kutta))]
    z = [runge_kutta[i][3] for i in range(len(runge_kutta))]
    ax.plot(x, y, z,
            label='parametric curve',
            color='red',
            marker='*',
            markersize=3)
    plt.legend(['Метод Эйлера', 'Метод Эйлера-Коши', 'Метод Рунге-Кутты'],
               loc='best')
    plt.title('График решения системы ОДУ')
    plt.show()


def equation_1(euler, euler_cauchy, runge_kutta):
    """
        Функция для построения графиков для 1 ур
        :param euler: list -- список с решением уравнения методом Эйлера
        вида: [x,y]
        :param euler_cauchy: list -- список с решением уравнения
        методом Эйлера-Коши вида: [x,y]
        :param runge_kutta: list -- список с решением уравнения
        методом Рунге-Кутты вида: [x,y]
        """
    x = [euler[i][1] for i in range(len(euler))]
    y = [euler[i][2] for i in range(len(euler))]
    plt.plot(x, y, color='pink', marker="*")
    x = [euler_cauchy[i][1] for i in range(len(euler_cauchy))]
    y = [euler_cauchy[i][2] for i in range(len(euler_cauchy))]
    plt.plot(x, y, marker='X',
             markersize=3,
             color='teal',
             markerfacecolor='greenyellow',
             markeredgecolor='teal')
    x = [runge_kutta[i][1] for i in range(len(runge_kutta))]
    y = [runge_kutta[i][2] for i in range(len(runge_kutta))]
    plt.plot(x, y, linestyle=':', linewidth=3, color='red')
    plt.title('График решения ОДУ', fontsize=18, fontname='Times New Roman')
    plt.legend(['Метод Эйлера',
                'Метод Эйлера-Коши',
                'Метод Рунге-Кутты'],
               loc='best')
    plt.show()


def enter():
    """
    Функция для ввода уравнений
    :return variable_new: list -- список переменных уравнений
    :return x0:float -- начальное условие x
    :return start: list -- начальные условия переменных в точке x0
    :return ab: list -- промежуток x
    :return f_eq: function _lambdifygenerated --
    :return variable_new1: list -- список переменных уравнений
    (переменные в формате y', z', y''...)
    необходим для ввода начальных условий при замене y'=y1
    :return N: количество уравнений в системе
    """
    N = int(input('Введите количество уравнений в системе: '))
    eq = []
    eq_new = []
    variable = []
    for i in range(N):
        eq_lot = input('Введите уравнение %d: ' % (i + 1))
        variable.append(eq_lot[0])
        eq.append(eq_lot)
    variable_new = []
    variable_new1 = []
    for i in range(N):
        eq_lot = eq[i]
        variable_new.append(variable[i])
        variable_new1.append(variable[i])
        f2s = []
        f1s = []
        for g in range(len(eq_lot)):
            if eq_lot[g] == '=':
                com = eq_lot[1:g]
                if (len(com) > 1):
                    for j in range(1, len(com)):
                        f1 = variable[i] + com[0] * j
                        f2 = variable[i] + str(j)
                        f2s.append(f2)
                        f1s.append(f1)
                        eq_new.append(f2)
                        variable_new.append(f2)
                        variable_new1.append(f1)
                    for r in range(N):
                        for k in range(len(f1s)):
                            eq[r] = re.sub(f1s[len(f1s) - k - 1],
                                           f2s[len(f2s) - k - 1], eq[r])
                eq_lot = eq[i]
                k = 0
                while (eq_lot[k] != '='):
                    k += 1
                eq_new.append(eq[i][k + 1:len(eq[i])])
                break
    x = symbols("x")
    f_eq = [0] * N
    f_eq = [lambdify([x, variable_new],
                     eq_new[i])for i in range(len(variable_new))]
    print('Введите начальные условия: ')
    x0 = float(input('x0= '))
    start = []
    for i in range(len(variable_new)):
        a = float(input('Значение %s в точке x0 ' % variable_new1[i]))
        start.append(a)
    b = float(input('Введите конец исследуемого промежутка: '))
    ab = [x0, b]
    print('Исследуемый промежуток:', ab)
    return variable_new, x0, start, ab, f_eq, variable_new1, N

