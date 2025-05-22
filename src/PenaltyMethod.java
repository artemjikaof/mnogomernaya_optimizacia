import java.util.Scanner;

public class PenaltyMethod {

    // Точное решение
    private static final double[] TRUE_X = {1.0 / 3, 1.0 / 3, 1.0 / 3};

    // Целевая функция
    static double f(double[] x) {
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2]; // x1^2 + x2^2 + x3^2
    }

    // Ограничение-равенство: g1(x) = x1 + x2 + x3 - 1
    static double g1(double[] x) {
        return x[0] + x[1] + x[2] - 1;
    }

    // Ограничение-неравенство: g2(x) = x1 - x2 - x3
    static double g2(double[] x) {
        return x[0] - x[1] - x[2];
    }

    // Штрафная функция
    static double F(double[] x, double r) {
        double penalty = 0.0;
        penalty += Math.pow(g1(x), 2);
        penalty += Math.pow(Math.max(0, g2(x)), 2);
        return f(x) + (r / 2.0) * penalty;
    }

    // Градиент штрафной функции
    static double[] gradF(double[] x, double r) {
        double[] grad = new double[3];

        grad[0] = 2*x[0];
        grad[1] = 2*x[1];
        grad[2] = 2*x[2];

        double dg1dx0 = 1;
        double dg1dx1 = 1;
        double dg1dx2 = 1;

        double g2val = g2(x);
        double dg2dx0 = 1;
        double dg2dx1 = -1;
        double dg2dx2 = -1;

        if (g2val > 0) {
            grad[0] += r * (dg1dx0 * g1(x) + dg2dx0 * g2val);
            grad[1] += r * (dg1dx1 * g1(x) + dg2dx1 * g2val);
            grad[2] += r * (dg1dx2 * g1(x) + dg2dx2 * g2val);
        } else {
            grad[0] += r * dg1dx0 * g1(x);
            grad[1] += r * dg1dx1 * g1(x);
            grad[2] += r * dg1dx2 * g1(x);
        }

        return grad;
    }

    // Матрица Гессе для штрафной функции
    static double[][] hessF(double[] x, double r) {
        double[][] H = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                H[i][j] = (i == j) ? 2 : 0;
            }
        }

        double dg1dx0 = 1;
        double dg1dx1 = 1;
        double dg1dx2 = 1;

        double g2val = g2(x);
        double dg2dx0 = 1;
        double dg2dx1 = -1;
        double dg2dx2 = -1;

        if (g2val > 0) {
            H[0][0] += r * (dg1dx0*dg1dx0 + dg2dx0*dg2dx0);
            H[0][1] += r * (dg1dx0*dg1dx1 + dg2dx0*dg2dx1);
            H[0][2] += r * (dg1dx0*dg1dx2 + dg2dx0*dg2dx2);

            H[1][0] += r * (dg1dx1*dg1dx0 + dg2dx1*dg2dx0);
            H[1][1] += r * (dg1dx1*dg1dx1 + dg2dx1*dg2dx1);
            H[1][2] += r * (dg1dx1*dg1dx2 + dg2dx1*dg2dx2);

            H[2][0] += r * (dg1dx2*dg1dx0 + dg2dx2*dg2dx0);
            H[2][1] += r * (dg1dx2*dg1dx1 + dg2dx2*dg2dx1);
            H[2][2] += r * (dg1dx2*dg1dx2 + dg2dx2*dg2dx2);
        } else {
            H[0][0] += r * dg1dx0*dg1dx0;
            H[0][1] += r * dg1dx0*dg1dx1;
            H[0][2] += r * dg1dx0*dg1dx2;

            H[1][0] += r * dg1dx1*dg1dx0;
            H[1][1] += r * dg1dx1*dg1dx1;
            H[1][2] += r * dg1dx1*dg1dx2;

            H[2][0] += r * dg1dx2*dg1dx0;
            H[2][1] += r * dg1dx2*dg1dx1;
            H[2][2] += r * dg1dx2*dg1dx2;
        }

        return H;
    }

    // Метод Ньютона
    static double[] newtonMinimize(double[] x0, double r, double eps1, double eps2, int M) {
        double[] x = x0.clone();
        int k = 0;

        while (k < M) {
            double[] grad = gradF(x, r);
            double normGrad = Math.sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
            if (normGrad <= eps1) break;

            double[][] H = hessF(x, r);
            double[][] Hinv = invertMatrix(H);
            if (Hinv == null) {
                x[0] -= 0.1 * grad[0];
                x[1] -= 0.1 * grad[1];
                x[2] -= 0.1 * grad[2];
                k++;
                continue;
            }

            double[] d = new double[3];
            for (int i = 0; i < 3; i++) {
                d[i] = -(Hinv[i][0] * grad[0] + Hinv[i][1] * grad[1] + Hinv[i][2] * grad[2]);
            }

            double[] xNew = new double[3];
            xNew[0] = x[0] + d[0];
            xNew[1] = x[1] + d[1];
            xNew[2] = x[2] + d[2];

            double dxNorm = Math.sqrt((xNew[0]-x[0])*(xNew[0]-x[0]) +
                    (xNew[1]-x[1])*(xNew[1]-x[1]) +
                    (xNew[2]-x[2])*(xNew[2]-x[2]));
            double df = Math.abs(F(xNew, r) - F(x, r));
            if (dxNorm < eps2 && df < eps2) break;

            x = xNew;
            k++;
        }
        return x;
    }

    // Вычисление обратной матрицы 3x3
    static double[][] invertMatrix(double[][] matrix) {
        double det = determinant(matrix);
        if (Math.abs(det) < 1e-10) return null;

        double[][] inv = new double[3][3];
        inv[0][0] =  (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / det;
        inv[0][1] = -(matrix[0][1] * matrix[2][2] - matrix[0][2] * matrix[2][1]) / det;
        inv[0][2] =  (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det;
        inv[1][0] = -(matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) / det;
        inv[1][1] =  (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det;
        inv[1][2] = -(matrix[0][0] * matrix[1][2] - matrix[0][2] * matrix[1][0]) / det;
        inv[2][0] =  (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / det;
        inv[2][1] = -(matrix[0][0] * matrix[2][1] - matrix[0][1] * matrix[2][0]) / det;
        inv[2][2] =  (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / det;
        return inv;
    }

    // Определитель матрицы 3x3
    static double determinant(double[][] m) {
        return m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
                - m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
                + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
    }

    // Расчёт абсолютной и относительной погрешности
    static void calculateErrors(double[] xApprox) {
        System.out.println("\nАнализ погрешности:");
        System.out.println("-------------------------------------------------------------------");
        System.out.printf("%6s | %10s | %10s | %10s | %12s\n", "Коорд.", "Точное", "Полученное", "Абс. погреш.", "Отн. погреш.");
        System.out.println("-------------------------------------------------------------------");

        for (int i = 0; i < 3; i++) {
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
        double[] x0 = new double[3];
        for (int i = 0; i < 3; i++) x0[i] = scanner.nextDouble();

        System.out.print("Введите r0: ");
        double r0 = scanner.nextDouble();

        System.out.print("Введите C: ");
        double C = scanner.nextDouble();

        System.out.print("Введите ε: ");
        double eps = scanner.nextDouble();

        double eps1 = 0.01;
        double eps2 = 0.01;
        int M = 100;

        int k = 0;
        double r = r0;
        double[] x = x0.clone();

        while (true) {
            x = newtonMinimize(x, r, eps1, eps2, M);

            double constraint1 = Math.abs(g1(x));
            double constraint2 = Math.max(0, g2(x));

            if (constraint1 <= eps && constraint2 <= eps) {
                break;
            }

            r *= C;
            k++;
        }

        System.out.printf("Точка минимума: (%.6f, %.6f, %.6f)\n", x[0], x[1], x[2]);
        System.out.printf("Значение функции в этой точке: %.6f\n", f(x));

        calculateErrors(x); // <-- Вызов метода расчёта погрешности
    }
}