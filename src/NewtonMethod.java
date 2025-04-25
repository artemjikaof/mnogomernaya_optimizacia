import java.util.Scanner;

public class NewtonMethod {

    // Функция f(x)
    public static double f(double[] x) {
        return Math.pow(x[0], 2) + Math.pow(x[1], 2) + x[0] * x[1] - 2 * x[0] - 6 * x[1];
    }

    // Градиент f(x)
    public static double[] gradient(double[] x) {
        double[] grad = new double[2];
        grad[0] = 2 * x[0] + x[1] - 2;
        grad[1] = 2 * x[1] + x[0] - 6;
        return grad;
    }

    // Матрица Гессе (в данном случае она постоянная)
    public static double[][] hessian() {
        return new double[][] { { 2, 1 }, { 1, 2 } };
    }

    // Обратная матрица Гессе
    public static double[][] inverseHessian() {
        double det = 3; // Определитель матрицы Гессе
        return new double[][] { { 2.0 / det, -1.0 / det }, { -1.0 / det, 2.0 / det } };
    }

    // Умножение матрицы на вектор
    public static double[] matrixVectorMultiply(double[][] matrix, double[] vector) {
        double[] result = new double[vector.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < vector.length; j++) {
                result[i] += matrix[i][j] * vector[j];
            }
        }
        return result;
    }

    // Вычисление нормы вектора
    public static double norm(double[] vector) {
        double sum = 0;
        for (double val : vector) {
            sum += val * val;
        }
        return Math.sqrt(sum);
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        // Ввод начальной точки
        System.out.println("Введите начальную точку x0 (x1, x2):");
        double x1 = scanner.nextDouble();
        double x2 = scanner.nextDouble();
        double[] x = { x1, x2 };

        // Ввод эпсилон-критериев
        System.out.println("Введите ε1:");
        double epsilon1 = scanner.nextDouble();
        System.out.println("Введите ε2:");
        double epsilon2 = scanner.nextDouble();

        // Ввод максимального числа итераций
        System.out.println("Введите предельное число итераций M:");
        int M = scanner.nextInt();

        // Матрица Гессе (постоянная)
        double[][] H = hessian();
        double[][] H_inv = inverseHessian();

        // Инициализация переменных
        int k = 0;
        double[] grad = gradient(x);
        double fx = f(x);
        double[] prevX = null;
        double prevFx = Double.NaN;

        // Главный цикл
        while (true) {
            // Проверка критерия окончания по градиенту
            if (norm(grad) <= epsilon1) {
                System.out.println("Расчет окончен: градиент меньше ε1.");
                break;
            }

            // Проверка предельного числа итераций
            if (k >= M) {
                System.out.println("Расчет окончен: достигнуто максимальное число итераций.");
                break;
            }

            // Вычисление направления поиска
            double[] d = matrixVectorMultiply(H_inv, grad); // d = -H_inv * grad
            for (int i = 0; i < d.length; i++) {
                d[i] = -d[i];
            }

            // Обновление точки
            double[] x_new = new double[2];
            x_new[0] = x[0] + d[0];
            x_new[1] = x[1] + d[1];

            // Вычисление нового значения функции
            double fx_new = f(x_new);

            // Проверка критериев остановки
            if (prevX != null && prevFx != Double.NaN) {
                if (norm(new double[] { x_new[0] - x[0], x_new[1] - x[1] }) < epsilon2 &&
                        Math.abs(fx_new - fx) < epsilon2) {
                    System.out.println("Расчет окончен: выполнены критерии остановки.");
                    break;
                }
            }

            // Обновление переменных
            prevX = x.clone();
            prevFx = fx;
            x = x_new;
            fx = fx_new;
            grad = gradient(x);

            // Вывод текущей информации
            System.out.printf("Итерация %d: x = (%.4f, %.4f), f(x) = %.4f\n", k, x[0], x[1], fx);

            k++;
        }

        // Вывод результата
        System.out.println("Минимум найден в точке: (" + x[0] + ", " + x[1] + ")");
        System.out.println("Значение функции в этой точке: " + f(x));
        scanner.close();
    }
}