import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Scanner;

public class OptimizationMethods1 {

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        // Ввод начальной точки x0
        System.out.print("Введите начальную точку x0 (через пробел): ");
        double x0_1 = scanner.nextDouble();
        double x0_2 = scanner.nextDouble();
        ArrayList<Double> x0 = new ArrayList<>(Arrays.asList(x0_1, x0_2));

        // Ввод параметров epsilon1, epsilon2 и maxIter
        System.out.print("Введите epsilon1: ");
        double epsilon1 = scanner.nextDouble();

        System.out.print("Введите epsilon2: ");
        double epsilon2 = scanner.nextDouble();

        System.out.print("Введите максимальное число итераций: ");
        int maxIter = scanner.nextInt();

        // Точное решение
        ArrayList<Double> exactSolution = new ArrayList<>(Arrays.asList(-2.0 / 3.0, 10.0 / 3.0));

        // Метод наискорейшего спуска
        ArrayList<Double> resultSteepest = fastestFall(x0, epsilon1, epsilon2, maxIter);
        printResult("Метод наискорейшего спуска", resultSteepest, exactSolution);

        // Метод сопряженных направлений
       // ArrayList<Double> resultConjugate = relatedDirections(x0, epsilon1, epsilon2, maxIter);
        //0printResult("Метод сопряженных направлений", resultConjugate, exactSolution);

        // Метод Ньютона
        ArrayList<Double> resultNewton = Newton(x0, epsilon1, epsilon2, maxIter);
        printResult("Метод Ньютона", resultNewton, exactSolution);

        scanner.close();
    }

    // Метод наискорейшего спуска
    public static ArrayList<Double> fastestFall(ArrayList<Double> x0, double e1, double e2, int maxIter) {
        ArrayList<Double> xSoZvedochkoy = new ArrayList<>();
        Map<Integer, ArrayList<Double>> x = new LinkedHashMap<>();
        x.put(0, x0);
        int k = 0;
        boolean flag = true;

        while (flag && k < maxIter) {
            if (norma1(gradFx(x.get(k).get(0), x.get(k).get(1))) < e1) {
                xSoZvedochkoy = x.get(k);
                flag = false;
            } else {
                double t = goldenRatioMethod(x.get(k), gradFx(x.get(k).get(0), x.get(k).get(1)));
                ArrayList<Double> xk = new ArrayList<>();
                xk.add(0, x.get(k).get(0) - t * gradFx(x.get(k).get(0), x.get(k).get(1)).get(0));
                xk.add(1, x.get(k).get(1) - t * gradFx(x.get(k).get(0), x.get(k).get(1)).get(1));
                x.put(k + 1, xk);

                if (norma(x.get(k + 1), x.get(k)) < e2 && Math.abs(f(x.get(k + 1)) - f(x.get(k))) < e2) {
                    xSoZvedochkoy = x.get(k + 1);
                    flag = false;
                } else {
                    k++;
                }
            }
        }

        return xSoZvedochkoy;
    }

    // Метод сопряженных направлений
    public static ArrayList<Double> relatedDirections(ArrayList<Double> x0, double e1, double e2, int maxIter) {
        ArrayList<Double> xSoZvedochkoy = new ArrayList<>();
        ArrayList<Double> x = new ArrayList<>(x0);
        Map<Integer, ArrayList<Double>> xi = new LinkedHashMap<>();
        xi.put(0, x);
        ArrayList<Double> d1 = new ArrayList<>(Arrays.asList(1.0, 0.0));
        ArrayList<Double> d2 = new ArrayList<>(Arrays.asList(0.0, 1.0));
        ArrayList<Double> d0 = d2;
        Map<Integer, ArrayList<Double>> d = new LinkedHashMap<>();
        d.put(1, d1);
        d.put(2, d2);
        d.put(0, d0);
        Map<Integer, ArrayList<Double>> y = new LinkedHashMap<>();
        y.put(0, xi.get(0));
        int n = 2;
        boolean flag = true;
        int i = 0;

        while (flag && i < maxIter) {
            int k = 0;
            double t = goldenRatioMethod(y.get(i), d.get(i + 1));
            ArrayList<Double> yi = new ArrayList<>();
            yi.add(0, y.get(i).get(0) + t * d.get(i + 1).get(0));
            yi.add(1, y.get(i).get(1) + t * d.get(i + 1).get(1));
            y.put(i + 1, yi);

            if (i < n - 1) {
                i++;
            } else {
                if (i == n - 1) {
                    if (norma(y.get(0), y.get(i + 1)) < e2) {
                        xSoZvedochkoy = y.get(i + 1);
                        flag = false;
                    } else {
                        xi.put(k + 1, y.get(i + 1));
                        if (norma(xi.get(k + 1), xi.get(k)) < e2) {
                            xSoZvedochkoy = xi.get(k + 1);
                            flag = false;
                        } else {
                            ArrayList<Double> dt = new ArrayList<>();
                            dt.add(0, y.get(i + 1).get(0) - y.get(1).get(0));
                            dt.add(1, y.get(i + 1).get(1) - y.get(1).get(1));
                            d.put(1, dt);
                            d.put(2, dt);
                            d.put(0, d.get(2));
                            i = 0;
                            y.put(0, xi.get(k + 1));
                            k++;
                        }
                    }
                }
            }
        }

        return xSoZvedochkoy;
    }

    // Метод Ньютона
    public static ArrayList<Double> Newton(ArrayList<Double> x0, double e1, double e2, int maxIter) {
        ArrayList<Double> xSoZvedochkoy = new ArrayList<>();
        Map<Integer, ArrayList<Double>> x = new LinkedHashMap<>();
        x.put(0, x0);
        int k = 0;
        boolean flag = true;
        ArrayList<Double> d = new ArrayList<>();

        while (flag && k < maxIter) {
            if (norma1(gradFx(x.get(k).get(0), x.get(k).get(1))) < e1) {
                xSoZvedochkoy = x.get(k);
                flag = false;
            } else {
                double t = 0;
                Map<Integer, ArrayList<Double>> gesseMatrix = GesseMatrix(x.get(k).get(0), x.get(k).get(1));
                Map<Integer, ArrayList<Double>> gesseRevMatrix = GesseRevMatrix(x.get(k).get(0), x.get(k).get(1));
                if (checkCharPlus(gesseRevMatrix)) {
                    d = multiplyMatrix(gradFx(x.get(k).get(0), x.get(k).get(1)), gesseRevMatrix);
                    t = 1;
                } else {
                    d.set(0, -gradFx(x.get(k).get(0), x.get(k).get(1)).get(0));
                    d.set(1, -gradFx(x.get(k).get(0), x.get(k).get(1)).get(1));
                    t = goldenRatioMethod(x.get(k), gradFx(x.get(k).get(0), x.get(k).get(1)));
                }

                ArrayList<Double> xk = new ArrayList<>();
                xk.add(x.get(k).get(0) + t * d.get(0));
                xk.add(x.get(k).get(1) + t * d.get(1));
                x.put(k + 1, xk);

                if (norma(x.get(k + 1), x.get(k)) < e2 && Math.abs(f(x.get(k + 1)) - f(x.get(k))) < e2) {
                    xSoZvedochkoy = x.get(k + 1);
                    flag = false;
                } else {
                    k++;
                }
            }
        }

        return xSoZvedochkoy;
    }

    // Функция
    private static double f(ArrayList<Double> x) {
        return x.get(0) * x.get(0) + x.get(1) * x.get(1) + x.get(0) * x.get(1) - 2 * x.get(0) - 6 * x.get(1);
    }

    // Градиент функции
    private static ArrayList<Double> gradFx(double x1, double x2) {
        ArrayList<Double> grad = new ArrayList<>();
        grad.add(2 * x1 + x2 - 2);
        grad.add(2 * x2 + x1 - 6);
        return grad;
    }

    // Гессиан функции
    private static Map<Integer, ArrayList<Double>> GesseMatrix(double x1, double x2) {
        Map<Integer, ArrayList<Double>> gesse = new LinkedHashMap<>();
        gesse.put(0, new ArrayList<>(Arrays.asList(2.0, 1.0)));
        gesse.put(1, new ArrayList<>(Arrays.asList(1.0, 2.0)));
        return gesse;
    }

    // Обратный Гессиан функции
    private static Map<Integer, ArrayList<Double>> GesseRevMatrix(double x1, double x2) {
        Map<Integer, ArrayList<Double>> gesseRev = new LinkedHashMap<>();
        double det = 2.0 * 2.0 - 1.0 * 1.0;
        gesseRev.put(0, new ArrayList<>(Arrays.asList(2.0 / det, -1.0 / det)));
        gesseRev.put(1, new ArrayList<>(Arrays.asList(-1.0 / det, 2.0 / det)));
        return gesseRev;
    }

    // Проверка положительной определенности матрицы
    private static boolean checkCharPlus(Map<Integer, ArrayList<Double>> matrix) {
        double det = matrix.get(0).get(0) * matrix.get(1).get(1) - matrix.get(0).get(1) * matrix.get(1).get(0);
        return det > 0;
    }

    // Умножение матрицы на вектор
    private static ArrayList<Double> multiplyMatrix(ArrayList<Double> vector, Map<Integer, ArrayList<Double>> matrix) {
        ArrayList<Double> result = new ArrayList<>();
        result.add(vector.get(0) * matrix.get(0).get(0) + vector.get(1) * matrix.get(0).get(1));
        result.add(vector.get(0) * matrix.get(1).get(0) + vector.get(1) * matrix.get(1).get(1));
        return result;
    }

    // Метод золотого сечения для нахождения оптимального шага t
    private static double goldenRatioMethod(ArrayList<Double> x, ArrayList<Double> grad) {
        double a = 0.0;
        double b = 1.0;
        double phi = (1 + Math.sqrt(5)) / 2;
        double eps = 0.01; // Точность для метода золотого сечения

        while (b - a > eps) {
            double x1 = b - (b - a) / phi;
            double x2 = a + (b - a) / phi;
            if (f(new ArrayList<>(Arrays.asList(x.get(0) - x1 * grad.get(0), x.get(1) - x1 * grad.get(1)))) <
                    f(new ArrayList<>(Arrays.asList(x.get(0) - x2 * grad.get(0), x.get(1) - x2 * grad.get(1))))) {
                b = x2;
            } else {
                a = x1;
            }
        }

        return (a + b) / 2;
    }

    // Норма вектора (евклидова)
    private static double norma(ArrayList<Double> x1, ArrayList<Double> x2) {
        return Math.sqrt(Math.pow(x1.get(0) - x2.get(0), 2) + Math.pow(x1.get(1) - x2.get(1), 2));
    }

    // Норма вектора (модуль)
    private static double norma1(ArrayList<Double> grad) {
        return Math.sqrt(grad.get(0) * grad.get(0) + grad.get(1) * grad.get(1));
    }

    // Вывод результата и погрешностей
    private static void printResult(String methodName, ArrayList<Double> result, ArrayList<Double> exactSolution) {
        double fx = f(result);
        double absError = norma(result, exactSolution);
        double relError = absError / norma(exactSolution, new ArrayList<>(Arrays.asList(0.0, 0.0)));

        System.out.printf("%s:\n", methodName);
        System.out.printf("x = (%.4f, %.4f)\n", result.get(0), result.get(1));
        System.out.printf("f(x) = %.4f\n", fx);
        System.out.printf("Абсолютная погрешность: %.4f\n", absError);
        System.out.printf("Относительная погрешность: %.4f\n", relError);
        System.out.println();
    }
}