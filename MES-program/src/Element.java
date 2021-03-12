public class Element {
    public int ID[];
    public double [][] Hl = new double[4][4];
    public double [][] Cl = new double[4][4];
    public double [] PC = new double[4];
    public void addtoHL(double [][]HBC){
        for (int n=0;n<4;n++){
            for(int m=0;m<4;m++){
                Hl[n][m]+=HBC[n][m];
            }
        }
    };
    public void setPC(double []PC){
      this.PC=PC;
    };
    public Element(){
        this.ID = new int[4];
    }
}

