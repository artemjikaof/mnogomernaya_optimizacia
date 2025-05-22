import java.util.Scanner;

public class PenaltyMethod2D {

    // Точное решение
    private static final double[] TRUE_X = {3.0, 1.0};

    // Целевая функция
    static double f(double[] x) {
        return 3 * x[0]*x[0] + 4 * x[0]*x[1] + 5 * x[1]*x[1]; // 3x1^2 + 4x1x2 + 5x2^2
    }

    // Ограничения (неравенства)
    static double g1(double[] x) { return -(x[0] + x[1] - 4); } // x1 + x2 >= 4 => -(x1 + x2 - 4) <= 0
    static double g2(double[] x) { return -x[0]; }             // x1 >= 0     => -x1 <= 0
    static double g3(double[] x) { return -x[1]; }             // x2 >= 0     => -x2 <= 0

    // Штрафная функция
    static double F(double[] x, double r) {
        double penalty = 0;
        penalty += Math.pow(Math.max(0, g1(x)), 2);
        penalty += Math.pow(Math.max(0, g2(x)), 2);
        penalty += Math.pow(Math.max(0, g3(x)), 2);
        return f(x) + (r / 2.0) * penalty;
    }

    // Градиент штрафной функции
    static double[] gradF(double[] x, double r) {
        double[] grad = new double[2];

        grad[0] = 6*x[0] + 4*x[1];
        grad[1] = 4*x[0] + 10*x[1];

        double dg1dx1 = -1;
        double dg1dx2 = -1;
        double dg2dx1 = -1;
        double dg2dx2 = 0;
        double dg3dx1 = 0;
        double dg3dx2 = -1;

        double g1val = g1(x);
        double g2val = g2(x);
        double g3val = g3(x);

        if (g1val > 0) {
            grad[0] += r * dg1dx1 * g1val;
            grad[1] += r * dg1dx2 * g1val;
        }
        if (g2val > 0) {
            grad[0] += r * dg2dx1 * g2val;
            grad[1] += r * dg2dx2 * g2val;
        }
        if (g3val > 0) {
            grad[0] += r * dg3dx1 * g3val;
            grad[1] += r * dg3dx2 * g3val;
        }

        return grad;
    }

    // Матрица Гессе для штрафной функции
    static double[][] hessF(double[] x, double r) {
        double[][] H = new double[2][2];

        H[0][0] = 6;
        H[0][1] = 4;
        H[1][0] = 4;
        H[1][1] = 10;

        double g1val = g1(x);
        double g2val = g2(x);
        double g3val = g3(x);

        if (g1val > 0) {
            H[0][0] += r * (-1)*(-1);
            H[0][1] += r * (-1)*(-1);
            H[1][0] += r * (-1)*(-1);
            H[1][1] += r * (-1)*(-1);
        }

        if (g2val > 0) {
            H[0][0] += r * (-1)*(-1);
        }

        if (g3val > 0) {
            H[1][1] += r * (-1)*(-1);
        }

        return H;
    }

    // Обратная матрица 2x2
    static double[][] invertMatrix2x2(double[][] matrix) {
        double det = matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
        if (Math.abs(det) < 1e-10) return null;

        double[][] inv = new double[2][2];
        inv[0][0] = matrix[1][1] / det;
        inv[0][1] = -matrix[0][1] / det;
        inv[1][0] = -matrix[1][0] / det;
        inv[1][1] = matrix[0][0] / det;

        return inv;
    }

    // Метод Ньютона
    static double[] newtonMinimize(double[] x0, double r, double eps1, double eps2, int M) {
        double[] x = x0.clone();
        int k = 0;

        while (k < M) {
            double[] grad = gradF(x, r);
            double normGrad = Math.sqrt(grad[0]*grad[0] + grad[1]*grad[1]);
            if (normGrad <= eps1) break;

            double[][] H = hessF(x, r);
            double[][] Hinv = invertMatrix2x2(H);
            if (Hinv == null) {
                x[0] -= 0.1 * grad[0];
                x[1] -= 0.1 * grad[1];
                k++;
                continue;
            }

            double[] d = new double[2];
            d[0] = -(Hinv[0][0] * grad[0] + Hinv[0][1] * grad[1]);
            d[1] = -(Hinv[1][0] * grad[0] + Hinv[1][1] * grad[1]);

            double[] xNew = new double[2];
            xNew[0] = x[0] + d[0];
            xNew[1] = x[1] + d[1];

            double dxNorm = Math.sqrt((xNew[0]-x[0])*(xNew[0]-x[0]) +
                    (xNew[1]-x[1])*(xNew[1]-x[1]));
            double df = Math.abs(F(xNew, r) - F(x, r));
            if (dxNorm < eps2 && df < eps2) break;

            x = xNew;
            k++;
        }
        return x;
    }

    // Расчёт абсолютной и относительной погрешности
    static void calculateErrors(double[] xApprox) {
        System.out.println("\nАнализ погрешности:");
        System.out.println("-------------------------------------------------------------------");
        System.out.printf("%6s | %10s | %10s | %10s | %12s\n", "Коорд.", "Точное", "Полученное", "Абс. погреш.", "Отн. погреш.");
        System.out.println("-------------------------------------------------------------------");

        for (int i = 0; i < 2; i++) {
            double exact = TRUE_X[i];
            double approx = xApprox[i];
            double absError = Math.abs(approx - exact);
            double relError = Math.abs(approx - exact) / Math.abs(exact);

            System.out.printf("  x[%d] | %10.6f | %10.6f | %10.6f | %10.4f\n",
                    i, exact, approx, absError, relError);
        }
        System.out.println("-------------------------------------------------------------------");
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        System.out.println("Введите начальную точку x0 (через пробел):");
        double[] x0 = new double[2];
        for (int i = 0; i < 2; i++) x0[i] = scanner.nextDouble();

        System.out.print("Введите r0: ");
        double r0 = scanner.nextDouble();

        System.out.print("Введите C: ");
        double C = scanner.nextDouble();

        System.out.print("Введите ε: ");
        double eps = scanner.nextDouble();

        double eps1 = 0.1;
        double eps2 = 0.1;
        int M = 100;

        int k = 0;
        double r = r0;
        double[] x = x0.clone();

        while (true) {
            x = newtonMinimize(x, r, eps1, eps2, M);

            boolean constraint1 = Math.max(0, g1(x)) <= eps;
            boolean constraint2 = Math.max(0, g2(x)) <= eps;
            boolean constraint3 = Math.max(0, g3(x)) <= eps;

            if (constraint1 && constraint2 && constraint3) {
                break;
            }

            r *= C;
            k++;
        }

        System.out.printf("Точка минимума: (%.6f, %.6f)\n", x[0], x[1]);
        System.out.printf("Значение функции в этой точке: %.6f\n", f(x));

        calculateErrors(x); // <-- Вызов метода расчёта погрешности
    }
}