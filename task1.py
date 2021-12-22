import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.misc import derivative


def f(x):
    """
    Функция преобразовывает введенную строку в математическое выражение
    с помощью функции eval().

    :param x: str -- строка
    :return: Математическое выражение

    """
    return eval(func,
                {'x': x,
                 'cos': math.cos, 'sin': math.sin, 'tg': math.tan,
                 'tan': math.tan, 'exp': math.exp, 'e': math.exp,
                 'log2': math.log2,
                 'ln': math.log, 'sqrt': math.sqrt})


def deriv(l_interval, r_interval, prec):
    """
    Функция для дифференцирования выражения.

    """
    # todo: поправить название переменной _
    n_1 = 0  # количество точек (int)

    deriv_x, deriv_y = [], []  # Значение точек х и
    # значение дифференцирования (list)

    # todo: поправить название переменной _
    for _ in np.arange(l_interval, r_interval, prec):
        try:
            df = (f(_ + prec) - f(_)) / prec  # дифференцирование
            # точек вручную (float)
        except ValueError:
            continue

        deriv_x.append(_)
        deriv_y.append(df)

        n_1 += 1

        try:
            df2 = derivative(f, _)  # дифференцирование встроенной
            # функцией scipy (float)
        except ValueError:
            print('{} | {} | {} | {} | {}'.format(n_1, round(_, 4),
                                                  df, None, None))
            continue

        print('{} | {} | {} | {} | {}'.format(n_1, round(_, 4), df,
                                              df2, df - df2))

    list_x, list_y = [], []  # список точек (list)

    # todo: поправить название переменной _
    for _ in np.arange(l_interval, r_interval, 0.01):
        try:
            list_y.append(f(_))
        except ValueError:
            continue
        list_x.append(_)

    plt.title("Графики функции и производной f(x)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.scatter(list_x, list_y, color='b')
    plt.scatter(list_x, list_y, label="f(x)", color='b')
    plt.scatter(deriv_x, deriv_y, color='c')
    plt.scatter(deriv_x, deriv_y, label="f'(x)", color='c')
    plt.legend()

    plt.show()


def integ(l_interval, r_interval, count_n):
    """
    Функция для интегрирования выражения.

    """
    s = 0  # площадь интегрируемой фигуры (float)
    value_h = (r_interval - l_interval) / count_n  # шаг сетки (float)
    for _ in np.linspace(l_interval, r_interval, count_n + 1)[1:-1]:
        try:
            s += f(_)
            s += f(l_interval) / 2
            s += f(r_interval) / 2
        except ValueError:
            continue

    s *= value_h
    print('Значение интеграла функции на отрезке '
          '[{}; {}] равно: {}'.format(l_interval, r_interval, s))


trig_func_lst = ['cos', 'sin', 'tg', 'sec', 'cosec']

func = input('Введите функцию: ')
func = func.replace('^', '**')

operation = int(input('Выберите один из вариантов:\n'
                      '0 - дифференцировать функцию,\n'
                      '1 - интегрировать функцию. \n'))

for element in trig_func_lst:
    if func.lower().find(element) != -1:
        angle = int(input('Введите: \n'
                          '0 - если углы в радианах, \n'
                          '1 - если в градусах\n'))

l_int = float(input('Введите левую границу интервала: '))
r_int = float(input('Введите правую границу интервала: '))
if l_int > r_int:
    print('Введенные данные некорректны, пожалуйста, повторите ввод.\n')
    quit()

pr = 1
count = 1
if operation == 0:
    step_pr = int(input('Введите: \n'
                        '0, если вы хотите задать точность, \n'
                        '1, если шаг.\n'))
    if step_pr == 0:
        pr = float(input('Введите точность: '))
        deriv(l_int, r_int, pr)
    elif step_pr == 1:
        pr = float(input('Введите шаг: '))
        deriv(l_int, r_int, pr)
    else:
        print('Введенные данные некорректны, пожалуйста, повторите ввод.\n')
        quit()
elif operation == 1:
    count_step = int(input('Введите: \n'
                           '0 - если вы хотите задать количество отрезков, \n'
                           '1 - если шаг.\n'))
    if count_step == 0:
        count = int(input('Введите кол-во отрезков: '))
        integ(l_int, r_int, count)
    elif count_step == 1:
        count = round((r_int - l_int) /
                      float(input('Введите шаг: ')))
        integ(l_int, r_int, count)
    else:
        print('Введенные данные некорректны, пожалуйста, повторите ввод.\n')
        quit()
else:
    print('Введенные данные некорректны, пожалуйста, повторите ввод.\n')
    quit()
