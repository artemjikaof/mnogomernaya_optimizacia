import java.util.*;

public class ConjugateDirections {
    static final double EPS1 = 0.01; // для интерполяции
    static final double EPS2 = 0.01; // для интерполяции
    static final double[] EXACT_SOLUTION = {-2.0 / 3.0, 10.0 / 3.0};

    public static double f(double[] x) {
        return x[0]*x[0] + x[1]*x[1] + x[0]*x[1] - 2*x[0] - 6*x[1];
    }

    public static double[] vectorAdd(double[] a, double[] b, double t) {
        return new double[] {a[0] + t*b[0], a[1] + t*b[1]};
    }

    public static double[] vectorSub(double[] a, double[] b) {
        return new double[] {a[0] - b[0], a[1] - b[1]};
    }

    public static double norm(double[] v) {
        return Math.sqrt(v[0]*v[0] + v[1]*v[1]);
    }

    public static double quadraticInterpolation(double[] y, double[] d) {
        double delta = 0.1;
        double x1 = 0;
        double x2 = x1 + delta;
        double f1 = f(vectorAdd(y, d, x1));
        double f2 = f(vectorAdd(y, d, x2));

        double x3 = (f1 > f2) ? x1 + 2*delta : x1 - delta;
        double f3 = f(vectorAdd(y, d, x3));

        while (true) {
            double numerator = 0.5*((x2*x2 - x3*x3)*f1 + (x3*x3 - x1*x1)*f2 + (x1*x1 - x2*x2)*f3);
            double denominator = (x2 - x3)*f1 + (x3 - x1)*f2 + (x1 - x2)*f3;

            if (Math.abs(denominator) < 1e-12) {
                x1 = (f1 <= f2 && f1 <= f3) ? x1 : (f2 <= f3 ? x2 : x3);
                x2 = x1 + delta;
                f1 = f(vectorAdd(y, d, x1));
                f2 = f(vectorAdd(y, d, x2));
                x3 = (f1 > f2) ? x1 + 2*delta : x1 - delta;
                f3 = f(vectorAdd(y, d, x3));
                continue;
            }

            double xBar = numerator / denominator;
            double fBar = f(vectorAdd(y, d, xBar));
            double Fmin = Math.min(Math.min(f1, f2), f3);
            double xMin = (Fmin == f1) ? x1 : (Fmin == f2) ? x2 : x3;

            if (Math.abs((Fmin - fBar) / fBar) < EPS1 && Math.abs((xMin - xBar) / xBar) < EPS2) {
                return xBar;
            }

            if (xBar >= Math.min(Math.min(x1, x2), x3) && xBar <= Math.max(Math.max(x1, x2), x3)) {
                double[] xs = {x1, x2, x3, xBar};
                Arrays.sort(xs);
                x1 = xs[0];
                x2 = xs[1];
                x3 = xs[2];
                f1 = f(vectorAdd(y, d, x1));
                f2 = f(vectorAdd(y, d, x2));
                f3 = f(vectorAdd(y, d, x3));
            } else {
                x1 = xBar;
                f1 = fBar;
                x2 = x1 + delta;
                f2 = f(vectorAdd(y, d, x2));
                x3 = (f1 > f2) ? x1 + 2*delta : x1 - delta;
                f3 = f(vectorAdd(y, d, x3));
            }
        }
    }

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        System.out.println("Введите начальную точку x0 (через пробел): ");
        double[] x0 = {sc.nextDouble(), sc.nextDouble()};

        System.out.print("Введите ε (точность): ");
        double epsilon = sc.nextDouble();

        System.out.print("Введите maxIter: ");
        int maxIter = sc.nextInt();

        double[] y = Arrays.copyOf(x0, 2);
        double[][] d = {{1, 0}, {0, 1}};
        int i = 0, k = 0;

        while (k < maxIter) {
            double[] direction = Arrays.copyOf(d[i], 2);
            double t = quadraticInterpolation(y, direction);
            double[] yNext = vectorAdd(y, direction, t);

            if (i < 1) {
                i++;
                y = yNext;
                continue;
            }

            // Проверка успешности поиска
            double[] checkDir = d[0];
            double tCheck = quadraticInterpolation(yNext, checkDir);
            double[] yCheck = vectorAdd(yNext, checkDir, tCheck);

            if (Math.abs(yCheck[0] - y[0]) < 1e-8 && Math.abs(yCheck[1] - y[1]) < 1e-8) {
                break;
            }

            // Построение нового направления
            double[] newDir = vectorSub(yNext, y);
            d[0] = d[1];
            d[1] = newDir;
            y = yNext;
            i = 0;
            k++;

            if (norm(newDir) < epsilon) break;
        }

        System.out.printf("Решение: x = (%.4f, %.4f)%n", y[0], y[1]);
        System.out.printf("Значение функции: f(x) = %.6f%n", f(y));
        double[] exact = EXACT_SOLUTION;
        double absError = Math.sqrt(Math.pow(y[0] - exact[0], 2) + Math.pow(y[1] - exact[1], 2));
        double relError = absError / Math.sqrt(exact[0]*exact[0] + exact[1]*exact[1]);
        System.out.printf("Абсолютная погрешность: %.6f%n", absError);
        System.out.printf("Относительная погрешность: %.6f%n", relError);
    }
}
