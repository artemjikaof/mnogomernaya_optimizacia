import java.util.Scanner;
import java.util.Arrays;

public class OptimizationMethods {

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        // Ввод начальной точки x0
        System.out.print("Введите начальную точку x0 (через пробел): ");
        double[] x0 = new double[2];
        for (int i = 0; i < 2; i++) {
            x0[i] = scanner.nextDouble();
        }

        // Ввод параметров epsilon1 и epsilon2
        System.out.print("Введите epsilon1: ");
        double epsilon1 = scanner.nextDouble();

        System.out.print("Введите epsilon2: ");
        double epsilon2 = scanner.nextDouble();

        // Ввод максимального числа итераций
        System.out.print("Введите максимальное число итераций: ");
        int maxIter = scanner.nextInt();

        // Точное решение
        double[] exactSolution = {-2.0 / 3.0, 10.0 / 3.0};
        double exactFunctionValue = -28.0 / 3.0;

        // Метод наискорейшего спуска
        Result resultSteepest = steepestDescent(x0, epsilon1, epsilon2, maxIter);
        printResult("Метод наискорейшего спуска", resultSteepest, exactSolution, exactFunctionValue);

        // Метод сопряженных направлений
        //Result resultConjugate = conjugateDirections(x0, epsilon1, epsilon2, maxIter);
        //printResult("Метод сопряженных направлений", resultConjugate, exactSolution, exactFunctionValue);

        // Метод Ньютона
        Result resultNewton = newtonMethod(x0, epsilon1, epsilon2, maxIter);
        printResult("Метод Ньютона", resultNewton, exactSolution, exactFunctionValue);

        scanner.close();
    }

    // Класс для хранения результата метода оптимизации
    static class Result {
        double[] x;
        int iterations;

        Result(double[] x, int iterations) {
            this.x = x;
            this.iterations = iterations;
        }
    }

    // Метод наискорейшего спуска
    public static Result steepestDescent(double[] x0, double epsilon1, double epsilon2, int maxIter) {
        double[] x = x0.clone();
        double[] grad = new double[2];
        double[] dx = new double[2];
        double fx, fxPrev;
        int k = 0;

        while (k < maxIter) {
            // Вычисляем градиент в текущей точке
            grad[0] = 2 * x[0] + x[1] - 2;
            grad[1] = 2 * x[1] + x[0] - 6;

            // Проверяем критерий окончания
            if (Math.sqrt(grad[0] * grad[0] + grad[1] * grad[1]) < epsilon1) {
                break;
            }

            // Вычисляем шаг t по методу квадратичной интерполяции
            double t = quadraticInterpolation(0.0, 1.0, epsilon1, epsilon2, new int[]{0}, grad, x);

            // Обновляем точку
            dx[0] = -t * grad[0];
            dx[1] = -t * grad[1];
            fxPrev = f(x);
            x[0] += dx[0];
            x[1] += dx[1];
            fx = f(x);

            // Проверяем условие окончания
            if (Math.sqrt(dx[0] * dx[0] + dx[1] * dx[1]) < epsilon2 && Math.abs(fx - fxPrev) < epsilon2) {
                break;
            }

            k++; // Увеличиваем счетчик итераций
        }

        return new Result(x, k);
    }

    // Метод сопряженных направлений
    public static Result conjugateDirections(double[] x0, double epsilon1, double epsilon2, int maxIter) {
        double[] x = x0.clone();
        double[] grad = new double[2];
        double[] d = new double[2];
        double[] y = new double[2];
        int n = 2; // Размерность пространства
        int i = 0;
        int k = 0;

        // Инициализация направлений
        d[0] = 1.0;
        d[1] = 0.0;

        while (k < maxIter) {
            // Вычисляем градиент в текущей точке
            grad[0] = 2 * x[0] + x[1] - 2;
            grad[1] = 2 * x[1] + x[0] - 6;

            // Проверяем критерий окончания
            if (Math.sqrt(grad[0] * grad[0] + grad[1] * grad[1]) < epsilon1) {
                break;
            }

            // Сохраняем текущую точку
            y[0] = x[0];
            y[1] = x[1];

            // Находим шаг t_i
            double t = quadraticInterpolation(0.0, 1.0, epsilon1, epsilon2, new int[]{0}, grad, x);

            // Обновляем точку
            x[0] -= t * grad[0];
            x[1] -= t * grad[1];

            // Проверяем успешность поиска
            if (i == n - 1) {
                if (Arrays.equals(y, x)) {
                    break;
                } else {
                    i = 0;
                }
            } else {
                i++;
            }

            // Проверяем условие окончания
            if (Math.sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1])) < epsilon2) {
                break;
            }

            k++; // Увеличиваем счетчик итераций
        }

        return new Result(x, k);
    }

    // Метод Ньютона
    public static Result newtonMethod(double[] x0, double epsilon1, double epsilon2, int maxIter) {
        double[] x = x0.clone();
        double[] grad = new double[2];
        double[][] hessian = new double[2][2];
        double[] dx = new double[2];
        double fx, fxPrev;
        int k = 0;

        while (k < maxIter) {
            // Вычисляем градиент
            grad[0] = 2 * x[0] + x[1] - 2;
            grad[1] = 2 * x[1] + x[0] - 6;

            // Проверяем критерий окончания
            if (Math.sqrt(grad[0] * grad[0] + grad[1] * grad[1]) < epsilon1) {
                break;
            }

            // Гессиан (в данном случае постоянен)
            hessian[0][0] = 2.0;
            hessian[0][1] = 1.0;
            hessian[1][0] = 1.0;
            hessian[1][1] = 2.0;

            // Решаем систему линейных уравнений H * dx = -grad
            dx = solveLinearSystem(hessian, new double[]{-grad[0], -grad[1]});

            // Обновляем точку
            fxPrev = f(x);
            x[0] += dx[0];
            x[1] += dx[1];
            fx = f(x);

            // Проверяем условие окончания
            if (Math.sqrt(dx[0] * dx[0] + dx[1] * dx[1]) < epsilon2 && Math.abs(fx - fxPrev) < epsilon2) {
                break;
            }

            k++; // Увеличиваем счетчик итераций
        }

        return new Result(x, k);
    }

    // Функция
    private static double f(double[] x) {
        return x[0] * x[0] + x[1] * x[1] + x[0] * x[1] - 2 * x[0] - 6 * x[1];
    }

    // Квадратичная интерполяция
    private static double quadraticInterpolation(double x1, double dx, double epsilon1, double epsilon2, int[] q, double[] grad, double[] x) {
        double xBar = 0.0;
        double x2 = x1 + dx;
        double f1 = f(new double[]{x[0] - x1 * grad[0], x[1] - x1 * grad[1]});
        double f2 = f(new double[]{x[0] - x2 * grad[0], x[1] - x2 * grad[1]});
        double x3, f3;
        boolean flag = true;
        int k = 0;
        if (f1 > f2) {
            x3 = x1 + 2 * dx;
        } else {
            x3 = x1 - dx;
        }

        f3 = f(new double[]{x[0] - x3 * grad[0], x[1] - x3 * grad[1]});

        while (flag) {
            k++;
            double Fmin = Math.min(f1, Math.min(f2, f3));
            double xmin = 0;
            if (f1 == Fmin) {
                xmin = x1;
            } else if (f2 == Fmin) {
                xmin = x2;
            } else {
                xmin = x3;
            }
            if (((x2 - x3) * f1 + (x3 - x1) * f2 + (x1 - x2) * f3) == 0) {
                x1 = xmin;
            } else {
                xBar = 0.5 * ((((Math.pow(x2, 2) - Math.pow(x3, 2)) * f1) + ((Math.pow(x3, 2) - Math.pow(x1, 2)) * f2) +
                        (((Math.pow(x1, 2) - Math.pow(x2, 2)) * f3))) /
                        (((x2 - x3) * f1) + ((x3 - x1) * f2) + ((x1 - x2) * f3)));
                double fXBar = f(new double[]{x[0] - xBar * grad[0], x[1] - xBar * grad[1]});
                if (Math.abs((Fmin - fXBar) / fXBar) < epsilon1 && Math.abs((xmin - xBar) / xBar) < epsilon2) {
                    q[0] = k;
                    flag = false;
                } else {
                    if (xBar <= x3 && xBar >= x1) {
                        if (fXBar < Fmin) {
                            if (x2 < xBar) {
                                x1 = x2;
                                x2 = xBar;
                                f1 = f2;
                                f2 = fXBar;
                            } else {
                                x3 = x2;
                                x2 = xBar;
                                f3 = f2;
                                f2 = fXBar;
                            }
                        } else {
                            if (x2 < xmin) {
                                x1 = x2;
                                x2 = xmin;
                                f1 = f2;
                                f2 = Fmin;
                            } else {
                                x3 = x2;
                                x2 = xmin;
                                f3 = f2;
                                f2 = Fmin;
                            }
                        }
                    } else {
                        x1 = xBar;
                        x2 = x1 + dx;
                        if (f(new double[]{x[0] - x1 * grad[0], x[1] - x1 * grad[1]}) > f(new double[]{x[0] - x2 * grad[0], x[1] - x2 * grad[1]}))
                            x3 = x1 + 2 * dx;
                        else
                            x3 = x1 - dx;
                        f1 = f(new double[]{x[0] - x1 * grad[0], x[1] - x1 * grad[1]});
                        f2 = f(new double[]{x[0] - x2 * grad[0], x[1] - x2 * grad[1]});
                        f3 = f(new double[]{x[0] - x3 * grad[0], x[1] - x3 * grad[1]});
                    }
                }

            }
        }
        return xBar;
    }

    // Решение системы линейных уравнений 2x2
    private static double[] solveLinearSystem(double[][] A, double[] b) {
        double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        double x1 = (A[1][1] * b[0] - A[0][1] * b[1]) / det;
        double x2 = (A[0][0] * b[1] - A[1][0] * b[0]) / det;
        return new double[]{x1, x2};
    }

    // Метод для вывода результатов
    private static void printResult(String methodName, Result result, double[] exactSolution, double exactFunctionValue) {
        double[] x = result.x;
        double functionValue = f(x);

        // Абсолютная погрешность для точки x
        double absErrorX = Math.sqrt(Math.pow(x[0] - exactSolution[0], 2) + Math.pow(x[1] - exactSolution[1], 2));

        // Относительная погрешность для точки x
        double relErrorX = absErrorX / Math.sqrt(exactSolution[0] * exactSolution[0] + exactSolution[1] * exactSolution[1]);

        // Абсолютная погрешность для значения функции
        double absErrorF = Math.abs(functionValue - exactFunctionValue);

        // Относительная погрешность для значения функции
        double relErrorF = absErrorF / Math.abs(exactFunctionValue);

        // Выводим результаты с количеством итераций k + 1
        System.out.printf("%s: x = (%f, %f), f(x) = %f, итераций: %d%n", methodName, x[0], x[1], functionValue, result.iterations + 1);
        System.out.printf("Абсолютная погрешность точки x: %f%n", absErrorX);
        System.out.printf("Относительная погрешность точки x: %f%n", relErrorX);
        System.out.printf("Абсолютная погрешность значения функции: %f%n", absErrorF);
        System.out.printf("Относительная погрешность значения функции: %f%n", relErrorF);
        System.out.println();
    }
}