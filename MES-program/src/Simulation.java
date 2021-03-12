import java.util.Arrays;


public class Simulation {
    public double[] count (double [][] H, double [][]C, double []P,Global global){
        int iteration=global.time/ global.dt;

        double [] temperature = new double [H.length];
        Arrays.fill(temperature, global.tempstart);
        for (int i = 0; i < iteration; i++) {
            double[][] h=new double[H.length][H.length];
            for (int j = 0; j <h.length ; j++) {
                for (int k = 0; k < h.length; k++) {
                    h[j][k]=H[j][k]+C[j][k]/ global.dt;
                }
            }
            double [] p = new double[P.length];
            for (int j = 0; j <h.length ; j++) {
                for (int k = 0; k < h.length; k++) {
                    p[j]+=(C[j][k]/ global.dt)*temperature[k];
                }
                p[j]-=P[j];
            }
            double [] newtemperature = CountMatrix.gauss(h,p);
            temperature=newtemperature.clone();
            Arrays.sort(newtemperature);
            System.out.println(newtemperature[0]+" "+newtemperature[newtemperature.length-1]);
        }
    return temperature;
    }
}
