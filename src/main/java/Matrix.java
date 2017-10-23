import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.TextArea;
import javafx.scene.paint.Color;

public class Matrix {

    private double[][] a;
    private double[][] backup;
    private double[] x = null;
    private double[][] statistic = null;
    private double[] r = null;
    private int size;
    private final double EPS = 1e-3;
    private final int MAXITER = 40;
    private int iters;

    public Matrix() {
        size = 4;
        //a = new double[][]{{1, 1000, 0, 8, 3}, {10, 4, 4, 5, 2}, {10, 4, 4, 5, 2}, {2, 4, 6, 5, 3}};
        a = new double[][]{{1, 1000, 0, 8, 3}, {10, 4, 4, 5, 2}, {-2, 8, -5, 0,2}, {2, 4, 6, 5, 3}};
        backup = new double[size][size+1];
        for (int i=0;i<size;i++) {
            for (int j=0;j<=size;j++)
                backup[i][j] = a[i][j];
        }
//        a = new double[][]{{1, 1, 1, 1}, {2, 1, 1, 2}, {1, 0, 1, -1}};
    }

    private void loadFromBackUp() {
        for (int i=0;i<size;i++) {
            for (int j=0;j<=size;j++)
                a[i][j] = backup[i][j];
        }
    }

    public void print() {
        for (int i=0;i<size;i++) {
            for (int j=0; j<size; j++) {
                if (a[i][j] > 0 && j>0) System.out.print("+");
                if (a[i][j] == 0) continue;
                if ((a[i][j] - (int)a[i][j]) < EPS) {
                    System.out.print((int)a[i][j] + "x"+(j+1));
                }
                else {
                    System.out.print(a[i][j] + "x"+(j+1));
                }
            }
            System.out.println("=" + a[i][size]);
        }
    }

    private void findX() {
        x = new double[size];
        for (int i=size-1; i>=0; i--) {
            double tmp = 0;
            for (int j = size-1; j>i; j--)
                tmp += x[j]*a[i][j];
            x[i] = (a[i][size] - tmp) / a[i][i];
        }
    }

    public void solveGauss() {
        loadFromBackUp();
        for (int k=0; k<size-1;k++) {
            gauusStep(k);
        }
        findX();
        //printSolve();
    }

    private int getNumberOfMaxedRow(int col) {
        double max = Math.abs(a[col][col]);
        int r = col;
        for (int i=col+1; i<size; i++) {
            if (Math.abs(a[i][col]) > max) {
                max = Math.abs(a[i][col]);
                r = i;
            }
        }
        return (max != 0)?r:-1;
    }

    private void swapRows(int m, int n) {
        for (int i=0;i<=size;i++) {
            double t = a[m][i];
            a[m][i] = a[n][i];
            a[n][i] = t;
        }
    }

    private void gauusStep(int row) {
        for (int i=row+1; i<size; i++) {
            double r = a[i][row] / a[row][row];
            for (int j=row; j<size;j++) {
                a[i][j] -= a[row][j]*r;
            }
            a[i][size] -= a[row][size]*r;
        }
    }

    public void solveGaussWithMainElement() throws BadSLAR {
        loadFromBackUp();
        for (int k=0; k<size;k++) {
            int h = getNumberOfMaxedRow(k);
            if (h == -1)
                throw new BadSLAR();
            if (h != k)
                swapRows(k, h);
            gauusStep(k);
        }
        findX();
    }

    public void printInMatrixMode() {
        if (x != null) {
            for (int i=0; i<size; i++) {
                System.out.print("|");
                for (int j=0; j<size;j++) {
                    System.out.print(String.format("%9.3f ",backup[i][j]));
                }
                System.out.print("| ");
                System.out.print((i == (size-1)/2)?"*":" ");
                System.out.print(" | ");
                System.out.print(String.format("%9.3f ",x[i]));
                System.out.print(" | ");
                System.out.print((i == (size-1)/2)?"=":" ");
                System.out.print(" | ");
                System.out.print(String.format("%9.3f ",backup[i][size]));
                System.out.println(" | ");
            }
        } else {
            System.out.println("Solve first");
        }
    }

    public void printSolve() {
        if (x != null) {
            for (int i=0; i<size; i++) {
                System.out.println("x"+(i+1)+" = " + x[i]);
            }
        } else {
            System.out.println("Solve first");
        }
    }

    public void calcVectorR() {
        if (x != null) {
            r = new double[size];
            for (int i=0; i<size;i++) {
                double s = 0;
                for (int j=0; j<size;j++) {
                    s += a[i][j] * x[j];
                }
                r[i] = s - a[i][size];
            }
        }
    }

    public void printR() {
        if (x != null) {
            calcVectorR();
            for (int i=0; i<size; i++) {
                System.out.println("r"+(i+1)+" = " + r[i]);
            }
        } else {
            System.out.println("Solve first");
        }
    }

    public void Nekrasov() {
        loadFromBackUp();
        statistic = new double[MAXITER][size+1];
        int iteration = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if ((Math.abs(backup[i][i]) > conv1()) || (Math.abs (backup[j][j]) > conv2())) {
                    statistic[iteration][size] = g();
                    while ((statistic[iteration][size] < EPS) && (iteration < MAXITER)) {
                        x[i] =(backup[i][size]- sum1()-sum2())/a[i][i];
                        statistic[iteration][i] = x[i];
                    }
                    iteration++;
                    System.out.println(iteration);
                }
            }
        }
        iters = iteration;
    }

    public double sum1() {
        double s=0;
        for(int i =0; i<size; i++){
            for (int j=1; j<i-1; j++){
                s +=backup[i][j]*x[j];
            }
        }
        return s;
    }

    public double sum2() {
        double s=0;
        for(int i =0; i<size; i++){
            for (int j=i+1; j<size; j++) {
                s +=backup[i][j]*x[j];
            }
        }
        return s;
    }

    public void drawStatistic(GraphicsContext gc) {
        gc.setFill(Color.GRAY);
        gc.fillRect(0,0,300, 300);
        gc.setStroke(Color.YELLOW);
        gc.setLineWidth(2);
        for(int i=0; i<iters; i++) {
        }
    }

    private double conv1() {
        double s = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j<size; j++){
                if (j != i)
                    s += Math.abs(backup[i][j]);
            }
        }
        return s;
    }

    private double conv2() {
        double s = 0;
        for (int j = 0; j < size; j++) {
            for (int i = 0; i < size; i++) {
                if (i != j)
                    s += Math.abs(backup[i][j]);
            }
        }
        return s;
    }


    private double g() {
        double g = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (Math.abs(x[i]) <= 1) g = Math.abs(x[i] - x[i-1]);
                else if (Math.abs(x[i]) > 1) g = Math.abs(x[i] - x[i-1]) / Math.abs(x[i]);
            }
        }
        return g;
    }

    public void printStatistic(TextArea area) {
        area.clear();
        StringBuilder sb = new StringBuilder();
        for(int i=0; i<iters; i++) {
            sb.append("x1=" + String.format("%9.3f", statistic[i][0]) + "x2=" + String.format("%9.3f", statistic[i][2]) +
                      "x3=" + String.format("%9.3f", statistic[i][2]) + "x4=" + String.format("%9.3f", statistic[i][3]) +
                      "\n");
            area.setText(sb.toString());
        }

    }
}
