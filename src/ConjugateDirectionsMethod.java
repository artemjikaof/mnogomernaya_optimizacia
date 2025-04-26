import java.util.Scanner;
//метод сопряженных направлений с золотым сечением
public class ConjugateDirectionsMethod {

    // Константа золотого сечения
    static final double GOLDEN_RATIO = (3 - Math.sqrt(5)) / 2; // 0.38196
    static final double INTERVAL_START = -10; // Левая граница интервала поиска шага
    static final double INTERVAL_END = 10;    // Правая граница интервала поиска шага

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        // Ввод начальной точки x0
        System.out.println("Введите начальную точку x0:");
        double[] x0 = new double[2];
        System.out.print("x1 = ");
        x0[0] = scanner.nextDouble();
        System.out.print("x2 = ");
        x0[1] = scanner.nextDouble();

        // Ввод точности ε
        System.out.print("Введите точность ε > 0: ");
        double epsilon = scanner.nextDouble();
        if (epsilon <= 0) {
            System.out.println("Ошибка: ε должно быть положительным числом.");
            return;
        }

        int n = x0.length; // Размерность задачи (n=2)

        // Формирование начальной системы направлений d1, d2, ..., dn и d0 = dn
        double[][] directions = new double[n + 1][n];
        for (int i = 0; i < n; i++) {
            directions[i][i] = 1.0;
        }
        directions[n] = directions[0].clone(); // d0 = dn

        int i = 0;        // Индекс текущего направления
        int k = 0;        // Номер итерации
        double[] y = x0.clone(); // Текущая точка y^0
        double[] xPrev = x0.clone(); // Предыдущая точка x^k
        double[] xNext;   // Следующая точка x^(k+1)

        // Главный цикл метода сопряжённых направлений
        while (true) {
            double[] direction = directions[i]; // Текущее направление поиска
            double t = goldenSectionSearch(y, direction, INTERVAL_START, INTERVAL_END, epsilon); // Поиск шага ti методом золотого сечения
            double[] yNext = vectorSum(y, scalarMultiply(direction, t)); // y^(i+1) = y^i + t_i * d_i

            if (i < n - 1) { // Пока не пройдены все направления
                y = yNext;
                i++;
                continue;
            }

            if (i == n - 1) { // После прохода всех направлений
                if (areVectorsEqual(yNext, x0, epsilon)) { // Проверка условия успешности
                    printSolution(yNext); // Вывод результата
                    return;
                } else {
                    y = yNext;
                    i++;
                    continue;
                }
            }

            if (i == n) { // Завершающий этап проверки
                if (areVectorsEqual(yNext, vectorSum(x0, scalarMultiply(directions[0], t)), epsilon)) {
                    printSolution(yNext);
                    return;
                } else {
                    xNext = yNext;
                    if (vectorNorm(vectorDiff(xNext, xPrev)) < epsilon) { // Условие останова по норме
                        printSolution(xNext);
                        return;
                    }

                    // Построение нового направления d0 = y^(n+1) - y^1
                    double[] newDirection = vectorDiff(yNext, vectorSum(x0, scalarMultiply(directions[0], t)));

                    // Формирование новой системы направлений
                    double[][] newDirections = new double[n + 1][n];
                    newDirections[0] = newDirection;
                    for (int j = 1; j < n; j++) {
                        newDirections[j] = directions[j].clone();
                    }
                    newDirections[n] = newDirections[0].clone();

                    // Проверка линейной независимости новой системы направлений
                    if (rank(newDirections, n) == n) {
                        directions = newDirections;
                    }

                    // Переход к следующей итерации
                    xPrev = xNext;
                    y = xNext.clone();
                    i = 0;
                    k++;
                }
            }
        }
    }

    // Функция f(x) = x1^2 + x2^2 + x1*x2 - 2x1 - 6x2
    static double f(double[] x) {
        double x1 = x[0];
        double x2 = x[1];
        return x1 * x1 + x2 * x2 + x1 * x2 - 2 * x1 - 6 * x2;
    }

    // Метод золотого сечения для одномерной минимизации вдоль направления d
    static double goldenSectionSearch(double[] y, double[] d, double a0, double b0, double l) {
        double a = a0;
        double b = b0;
        double y0 = a + GOLDEN_RATIO * (b - a);
        double z0 = a + b - y0;

        double fy = f(vectorSum(y, scalarMultiply(d, y0)));
        double fz = f(vectorSum(y, scalarMultiply(d, z0)));

        while (Math.abs(b - a) > l) {
            if (fy <= fz) {
                b = z0;
                z0 = y0;
                fz = fy;
                y0 = a + b - y0;
                fy = f(vectorSum(y, scalarMultiply(d, y0)));
            } else {
                a = y0;
                y0 = z0;
                fy = fz;
                z0 = a + b - z0;
                fz = f(vectorSum(y, scalarMultiply(d, z0)));
            }
        }
        return (a + b) / 2.0; // Приближённое значение шага t
    }

    // Умножение вектора на скаляр
    static double[] scalarMultiply(double[] v, double scalar) {
        double[] res = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            res[i] = v[i] * scalar;
        }
        return res;
    }

    // Сумма двух векторов
    static double[] vectorSum(double[] a, double[] b) {
        double[] res = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            res[i] = a[i] + b[i];
        }
        return res;
    }

    // Разность двух векторов
    static double[] vectorDiff(double[] a, double[] b) {
        double[] res = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            res[i] = a[i] - b[i];
        }
        return res;
    }

    // Норма вектора (евклидова)
    static double vectorNorm(double[] v) {
        double sum = 0.0;
        for (double val : v) {
            sum += val * val;
        }
        return Math.sqrt(sum);
    }

    // Проверка равенства двух векторов с заданной точностью ε
    static boolean areVectorsEqual(double[] a, double[] b, double eps) {
        for (int i = 0; i < a.length; i++) {
            if (Math.abs(a[i] - b[i]) > eps) {
                return false;
            }
        }
        return true;
    }

    // Вывод оптимальной точки и значения функции в ней
    static void printSolution(double[] x) {
        System.out.printf("Точка минимума: (%.5f, %.5f)%n", x[0], x[1]);
        System.out.printf("Значение функции в точке минимума: %.5f%n", f(x));
    }

    // Вычисление ранга матрицы (метод Гаусса)
    static int rank(double[][] matrix, int n) {
        double[][] temp = new double[n][n];
        for (int i = 0; i < n; i++) {
            temp[i] = matrix[i].clone();
        }

        int r = 0;
        for (int i = 0; i < n; i++) {
            int pivot = -1;
            for (int j = i; j < n; j++) {
                if (Math.abs(temp[j][i]) > 1e-8) {
                    pivot = j;
                    break;
                }
            }
            if (pivot == -1) {
                continue;
            }

            r++;
            double[] tmp = temp[i];
            temp[i] = temp[pivot];
            temp[pivot] = tmp;

            for (int j = i + 1; j < n; j++) {
                double factor = temp[j][i] / temp[i][i];
                for (int k = i; k < n; k++) {
                    temp[j][k] -= factor * temp[i][k];
                }
            }
        }
        return r;
    }
}
