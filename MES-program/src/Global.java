import java.io.*;

import java.util.Scanner;

public class Global {
    int pc;
    double W; //dlugosc
    double H; //wysokosc
    int nH;
    int nW;
    int nE;
    int nN;
    double p;
    double cp;
    double otoczenie;
    double alfa;
    int tempstart;
    int time;
    int dt;


    public Global(){
        try{
            File file = new File("dane.txt"); // tu zmienic na dane2.txt jak dla drugiego test case
            Scanner Reader = new Scanner(file);
            this.pc= Reader.nextInt();
            this.W = Reader.nextDouble();
            this.H = Reader.nextDouble();
            this.nW = Reader.nextInt();
            this.nH = Reader.nextInt();
            this.p = Reader.nextDouble();
            this.cp= Reader.nextDouble();
            this.nE = (nW-1)*(nH-1);
            this.nN = nW*nH;
            this.otoczenie= Reader.nextDouble();
            this.alfa= Reader.nextDouble();
            this.tempstart=Reader.nextInt();
            this.time= Reader.nextInt();
            this.dt= Reader.nextInt();
            Reader.close();
        }
        catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
