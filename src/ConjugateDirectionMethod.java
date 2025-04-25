import java.util.Scanner;

public class ConjugateDirectionMethod {

    public static void main(String[] args) {
        testQuadraticInterpolation();

        Scanner scanner = new Scanner(System.in);

        // Ввод начальной точки
        System.out.println("Введите начальную точку x0 (через пробел):");
        double x0 = scanner.nextDouble();
        double y0 = scanner.nextDouble();

        // Ввод точности
        System.out.println("Введите точность:");
        double epsilon = scanner.nextDouble();

        // Начальные параметры
        double[] startPoint = {x0, y0};
        int dimension = startPoint.length;
        double[][] directions = new double[dimension][dimension];

        // Инициализация направлений
        for (int i = 0; i < dimension; i++) {
            directions[i][i] = 1.0;
        }

        double[] currentPoint = startPoint.clone();
        double[] previousPoint = startPoint.clone();
        double[] yi = startPoint.clone();
        double[] yNext = new double[dimension];
        double[] conjugateDirection = new double[dimension];
        double[] newDirections = new double[dimension];
        int iterationCount = 0;
        boolean successSearch = false;

        while (true) {
            successSearch = false;
            for (int i = 0; i < dimension; i++) {
                double ti = findStep(yi, directions[i], epsilon);
                for (int j = 0; j < dimension; j++) {
                    yNext[j] = yi[j] + ti * directions[i][j];
                }
                if (i < dimension - 1) {
                    System.arraycopy(yNext, 0, yi, 0, dimension);
                } else {
                    successSearch = checkSuccessSearch(yNext, yi, dimension);
                    if (successSearch) {
                        break;
                    }
                    System.arraycopy(yNext, 0, yi, 0, dimension);
                }
            }

            if (successSearch) {
                break;
            }

            // Построение сопряженного направления
            for (int i = 0; i < dimension; i++) {
                conjugateDirection[i] = yNext[i] - yi[i];
            }

            // Проверка линейной зависимости
            if (isLinearlyIndependent(conjugateDirection, directions)) {
                System.arraycopy(conjugateDirection, 0, newDirections, 0, dimension);
                System.arraycopy(directions, 0, newDirections, 1, dimension - 1);
                System.arraycopy(newDirections, 0, directions, 0, dimension);
            }

            // Обновление текущей точки
            System.arraycopy(yNext, 0, currentPoint, 0, dimension);

            // Проверка условия окончания
            if (norm(currentPoint, previousPoint) < epsilon) {
                break;
            }

            System.arraycopy(currentPoint, 0, previousPoint, 0, dimension);
            iterationCount++;
        }

        System.out.println("Оптимальная точка: (" + currentPoint[0] + ", " + currentPoint[1] + ")");
        System.out.println("Количество итераций: " + iterationCount);
    }

    private static double findStep(double[] point, double[] direction, double epsilon) {
        double x1 = 0.0;
        double dx = 0.1;
        int[] q = new int[1];
        return quadraticInterpolation(x1, dx, epsilon, epsilon, q, point, direction);
    }

    private static double quadraticInterpolation(double x1, double dx, double epsilon1, double epsilon2, int[] q, double[] point, double[] direction) {
        double xBar = 0.0;
        double x2 = x1 + dx;
        double f1 = f(point, direction, x1);
        double f2 = f(point, direction, x2);
        double x3, f3;
        boolean flag = true;
        int k = 0;
        if (f1 > f2) {
            x3 = x1 + 2 * dx;
        } else {
            x3 = x1 - dx;
        }

        f3 = f(point, direction, x3);

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
                double fXBar = f(point, direction, xBar);
                if (Math.abs((Fmin - fXBar) / fXBar) < epsilon1 && Math.abs((xmin - xBar) / xBar) < epsilon2) {
                    q[0] = k;
                    flag = false;
                } else {
                    if (xBar <= x3 & xBar >= x1) {
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
                        if (f(point, direction, x1) > f(point, direction, x2))
                            x3 = x1 + 2 * dx;
                        else
                            x3 = x1 - dx;
                        f1 = f(point, direction, x1);
                        f2 = f(point, direction, x2);
                        f3 = f(point, direction, x3);
                    }
                }
            }
        }
        return xBar;
    }

    private static double f(double[] point, double[] direction, double t) {
        // Используем тестовую функцию для проверки
        return fTest(point, direction, t);
    }

    private static double fTest(double[] point, double[] direction, double t) {
        double x = point[0] + t * direction[0];
        return (x - 3) * (x - 3);
    }

    private static boolean checkSuccessSearch(double[] yNext, double[] yi, int dimension) {
        for (int i = 0; i < dimension; i++) {
            if (yNext[i] != yi[i]) {
                return false;
            }
        }
        return true;
    }

    private static boolean isLinearlyIndependent(double[] newDirection, double[][] directions) {
        double[][] matrix = new double[directions.length + 1][directions.length + 1];
        for (int i = 0; i < directions.length; i++) {
            System.arraycopy(directions[i], 0, matrix[i], 0, directions.length);
        }
        System.arraycopy(newDirection, 0, matrix[matrix.length - 1], 0, matrix.length - 1);
        return determinant(matrix) != 0;
    }

    private static double determinant(double[][] matrix) {
        int n = matrix.length;
        if (n == 1) {
            return matrix[0][0];
        }
        double det = 0;
        for (int i = 0; i < n; i++) {
            det += Math.pow(-1, i) * matrix[0][i] * determinant(minor(matrix, 0, i));
        }
        return det;
    }

    private static double[][] minor(double[][] matrix, int row, int col) {
        int n = matrix.length;
        double[][] minor = new double[n - 1][n - 1];
        for (int i = 0, mi = 0; i < n; i++) {
            if (i == row) continue;
            for (int j = 0, mj = 0; j < n; j++) {
                if (j == col) continue;
                minor[mi][mj++] = matrix[i][j];
            }
            mi++;
        }
        return minor;
    }

    private static double norm(double[] x, double[] y) {
        double sum = 0;
        for (int i = 0; i < x.length; i++) {
            sum += Math.pow(x[i] - y[i], 2);
        }
        return Math.sqrt(sum);
    }

    public static void testQuadraticInterpolation() {
        double x1 = 0.0;
        double dx = 0.1;
        double epsilon1 = 1e-6;
        double epsilon2 = 1e-6;
        int[] q = new int[1];
        double[] point = {0, 0}; // Не важно для тестирования f(t)
        double[] direction = {1, 0}; // Направление по оси x

        double result = quadraticInterpolation(x1, dx, epsilon1, epsilon2, q, point, direction);
        System.out.println("Результат квадратичной интерполяции: " + result);
    }
}