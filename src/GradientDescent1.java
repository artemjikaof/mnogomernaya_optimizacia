import java.util.Scanner;
// метод наискорейшего спуска с методом золотого сечения
public class GradientDescent1 {

    // Константы по умолчанию
    private static final double GOLDEN_RATIO = (3 - Math.sqrt(5)) / 2; // 0.38196

    // Функция f(x)
    public static double function(double[] x) {
        return x[0] * x[0] + x[1] * x[1] + x[0] * x[1] - 2 * x[0] - 6 * x[1];
    }

    // Градиент функции f(x)
    public static double[] gradient(double[] x) {
        double[] grad = new double[2];
        grad[0] = 2 * x[0] + x[1] - 2;
        grad[1] = 2 * x[1] + x[0] - 6;
        return grad;
    }

    // Метод золотого сечения для поиска оптимального шага
    public static double goldenSectionSearch(double[] x, double[] grad, double a, double b, double l) {
        double y = a + GOLDEN_RATIO * (b - a);
        double z = a + b - y;
        double fy = function(new double[]{x[0] - y * grad[0], x[1] - y * grad[1]});
        double fz = function(new double[]{x[0] - z * grad[0], x[1] - z * grad[1]});

        while (Math.abs(b - a) > l) {
            if (fy <= fz) {
                b = z;
                z = y;
                y = a + b - z;
                fz = fy;
                fy = function(new double[]{x[0] - y * grad[0], x[1] - y * grad[1]});
            } else {
                a = y;
                y = z;
                z = a + b - y;
                fy = fz;
                fz = function(new double[]{x[0] - z * grad[0], x[1] - z * grad[1]});
            }
        }
        return (a + b) / 2;
    }

    // Основной метод оптимизации
    public static double[] optimize(double[] x, double epsilon1, double epsilon2, int maxIterations) {
        double[] grad;
        double tk;
        int k = 0;
        double normGrad;
        double normXDiff;
        double funcDiff;
        double prevFuncValue = function(x);

        do {
            grad = gradient(x);
            normGrad = Math.sqrt(grad[0] * grad[0] + grad[1] * grad[1]);

            if (normGrad < epsilon1) {
                break;
            }

            tk = goldenSectionSearch(x, grad, 0, 1, epsilon2);
            double[] newX = {x[0] - tk * grad[0], x[1] - tk * grad[1]};
            normXDiff = Math.sqrt((newX[0] - x[0]) * (newX[0] - x[0]) + (newX[1] - x[1]) * (newX[1] - x[1]));
            funcDiff = Math.abs(function(newX) - prevFuncValue);
            prevFuncValue = function(newX);

            if (normXDiff < epsilon2 && funcDiff < epsilon2) {
                x = newX;
                break;
            }

            x = newX;
            k++;
        } while (k < maxIterations);

        return x;
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        // Ввод начальной точки
        System.out.print("Введите начальную точку x0 (x1, x2):");
        double x1 = scanner.nextDouble();
        double x2 = scanner.nextDouble();
        double[] x = {x1, x2};

        // Ввод параметров epsilon1 и epsilon2
        System.out.print("Введите epsilon1: ");
        double epsilon1 = scanner.nextDouble();
        System.out.print("Введите epsilon2: ");
        double epsilon2 = scanner.nextDouble();

        // Ввод максимального числа итераций
        System.out.print("Введите максимальное число итераций: ");
        int maxIterations = scanner.nextInt();

        // Оптимизация
        double[] result = optimize(x, epsilon1, epsilon2, maxIterations);
        System.out.println("Точка минимума: (" + result[0] + ", " + result[1] + ")");
        System.out.println("Значение функции в точке минимума: " + function(result));

        scanner.close();
    }
}