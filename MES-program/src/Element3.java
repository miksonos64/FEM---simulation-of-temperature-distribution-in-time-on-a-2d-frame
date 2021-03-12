import java.util.List;


public class Element3 implements IElem3{

    double [] E ;
    double [] N ;
    double [] scale1 ;
    double [] scale2 ;
    int numberP;

    double L[][] ;


    double pktcalkE[][];

    double pktcalkN[][];

    public Element3 (int numberP) {
        this.numberP=numberP;
        if (numberP==2) {
            E =new double[] {-(1/Math.sqrt(3)), 1/Math.sqrt(3), 1/Math.sqrt(3), -(1/Math.sqrt(3))};
            N =new double[] {-(1/Math.sqrt(3)), -(1/Math.sqrt(3)), 1/Math.sqrt(3), 1/Math.sqrt(3)};
            scale1 =new double[] {1, 1, 1, 1};
            scale2 =new double[] {1, 1, 1, 1};


            L = new double[][]{
                    {0.25 * (1 - E[0]) * (1 - N[0]), 0.25 * (1 + E[0]) * (1 - N[0]), 0.25 * (1 + E[0]) * (1 + N[0]), 0.25 * (1 - E[0]) * (1 + N[0])},
                    {0.25 * (1 - E[1]) * (1 - N[1]), 0.25 * (1 + E[1]) * (1 - N[1]), 0.25 * (1 + E[1]) * (1 + N[1]), 0.25 * (1 - E[1]) * (1 + N[1])},
                    {0.25 * (1 - E[2]) * (1 - N[2]), 0.25 * (1 + E[2]) * (1 - N[2]), 0.25 * (1 + E[2]) * (1 + N[2]), 0.25 * (1 - E[2]) * (1 + N[2])},
                    {0.25 * (1 - E[3]) * (1 - N[3]), 0.25 * (1 + E[3]) * (1 - N[3]), 0.25 * (1 + E[3]) * (1 + N[3]), 0.25 * (1 - E[3]) * (1 + N[3])}
            };

            pktcalkE =new double[][]
                    {
                            {-0.25 * (1 - N[0]), 0.25 * (1 - N[0]), 0.25 * (1 + N[0]), -0.25 * (1 + N[0])},
                            {-0.25 * (1 - N[1]), 0.25 * (1 - N[1]), 0.25 * (1 + N[1]), -0.25 * (1 + N[1])},
                            {-0.25 * (1 - N[2]), 0.25 * (1 - N[2]), 0.25 * (1 + N[2]), -0.25 * (1 + N[2])},
                            {-0.25 * (1 - N[3]), 0.25 * (1 - N[3]), 0.25 * (1 + N[3]), -0.25 * (1 + N[3])}
                    };

            pktcalkN =new double[][]
                    {
                            {-0.25 * (1 - E[0]), -0.25 * (1 + E[0]), 0.25 * (1 + E[0]), 0.25 * (1 - E[0])},
                            {-0.25 * (1 - E[1]), -0.25 * (1 + E[1]), 0.25 * (1 + E[1]), 0.25 * (1 - E[1])},
                            {-0.25 * (1 - E[2]), -0.25 * (1 + E[2]), 0.25 * (1 + E[2]), 0.25 * (1 - E[2])},
                            {-0.25 * (1 - E[3]), -0.25 * (1 + E[3]), 0.25 * (1 + E[3]), 0.25 * (1 - E[3])}
                    };
        }
        else if (numberP==3){
            E = new double[]{-1*Math.sqrt(3.0/5.0), 0, Math.sqrt(3.0/5.0), -1*Math.sqrt(3.0/5.0), 0, Math.sqrt(3.0/5.0), -1*Math.sqrt(3.0/5.0), 0, Math.sqrt(3.0/5.0)};
            N = new double[]{-1*Math.sqrt(3.0/5.0), -1*Math.sqrt(3.0/5.0), -1*Math.sqrt(3.0/5.0), 0,0,0, Math.sqrt(3.0/5.0), Math.sqrt(3.0/5.0), Math.sqrt(3.0/5.0)};
            scale1 =new double[] {5.0/9.0, 8.0/9.0, 5.0/9.0, 5.0/9.0, 8.0/9.0, 5.0/9.0, 5.0/9.0, 8.0/9.0, 5.0/9.0};
            scale2 =new double[] {5.0/9.0, 5.0/9.0, 5.0/9.0, 8.0/9.0, 8.0/9.0, 8.0/9.0, 5.0/9.0, 5.0/9.0, 5.0/9.0};
            pktcalkE =new double[][]
                    {
                            { -0.25*(1 - N[0]) , 0.25 * (1 - N[0]) , 0.25 * (1 + N[0]) , -0.25 * (1 + N[0]) },
                            { -0.25*(1 - N[1]) , 0.25 * (1 - N[1]) , 0.25 * (1 + N[1]) , -0.25 * (1 + N[1]) },
                            { -0.25*(1 - N[2]) , 0.25 * (1 - N[2]) , 0.25 * (1 + N[2]) , -0.25 * (1 + N[2]) },
                            { -0.25*(1 - N[3]) , 0.25 * (1 - N[3]) , 0.25 * (1 + N[3]) , -0.25 * (1 + N[3]) },
                            { -0.25*(1 - N[4]) , 0.25 * (1 - N[4]) , 0.25 * (1 + N[4]) , -0.25 * (1 + N[4]) },
                            { -0.25*(1 - N[5]) , 0.25 * (1 - N[5]) , 0.25 * (1 + N[5]) , -0.25 * (1 + N[5]) },
                            { -0.25*(1 - N[6]) , 0.25 * (1 - N[6]) , 0.25 * (1 + N[6]) , -0.25 * (1 + N[6]) },
                            { -0.25*(1 - N[7]) , 0.25 * (1 - N[7]) , 0.25 * (1 + N[7]) , -0.25 * (1 + N[7]) },
                            { -0.25*(1 - N[8]) , 0.25 * (1 - N[8]) , 0.25 * (1 + N[8]) , -0.25 * (1 + N[8]) }
                    };

            pktcalkN =new double[][]
                    {
                            { -0.25*(1 - E[0]) , -0.25 * (1 + E[0]) , 0.25 * (1 + E[0]) , 0.25 * (1 - E[0]) },
                            { -0.25*(1 - E[1]) , -0.25 * (1 + E[1]) , 0.25 * (1 + E[1]) , 0.25 * (1 - E[1]) },
                            { -0.25*(1 - E[2]) , -0.25 * (1 + E[2]) , 0.25 * (1 + E[2]) , 0.25 * (1 - E[2]) },
                            { -0.25*(1 - E[3]) , -0.25 * (1 + E[3]) , 0.25 * (1 + E[3]) , 0.25 * (1 - E[3]) },
                            { -0.25*(1 - E[4]) , -0.25 * (1 + E[4]) , 0.25 * (1 + E[4]) , 0.25 * (1 - E[4]) },
                            { -0.25*(1 - E[5]) , -0.25 * (1 + E[5]) , 0.25 * (1 + E[5]) , 0.25 * (1 - E[5]) },
                            { -0.25*(1 - E[6]) , -0.25 * (1 + E[6]) , 0.25 * (1 + E[6]) , 0.25 * (1 - E[6]) },
                            { -0.25*(1 - E[7]) , -0.25 * (1 + E[7]) , 0.25 * (1 + E[7]) , 0.25 * (1 - E[7]) },
                            { -0.25*(1 - E[8]) , -0.25 * (1 + E[8]) , 0.25 * (1 + E[8]) , 0.25 * (1 - E[8]) }
                    };

             L =new double[][]{
                    {0.25 *(1-E[0]) * (1 - N[0]),0.25 *(1+E[0]) * (1 - N[0]),0.25 *(1+E[0]) * (1 + N[0]),0.25 *(1-E[0]) * (1 + N[0])},
                    {0.25 *(1-E[1]) * (1 - N[1]),0.25 *(1+E[1]) * (1 - N[1]),0.25 *(1+E[1]) * (1 + N[1]),0.25 *(1-E[1]) * (1 + N[1])},
                    {0.25 *(1-E[2]) * (1 - N[2]),0.25 *(1+E[2]) * (1 - N[2]),0.25 *(1+E[2]) * (1 + N[2]),0.25 *(1-E[2]) * (1 + N[2])},
                    {0.25 *(1-E[3]) * (1 - N[3]),0.25 *(1+E[3]) * (1 - N[3]),0.25 *(1+E[3]) * (1 + N[3]),0.25 *(1-E[3]) * (1 + N[3])},
                    {0.25 *(1-E[4]) * (1 - N[4]),0.25 *(1+E[4]) * (1 - N[4]),0.25 *(1+E[4]) * (1 + N[4]),0.25 *(1-E[4]) * (1 + N[4])},
                    {0.25 *(1-E[5]) * (1 - N[5]),0.25 *(1+E[5]) * (1 - N[5]),0.25 *(1+E[5]) * (1 + N[5]),0.25 *(1-E[5]) * (1 + N[5])},
                    {0.25 *(1-E[6]) * (1 - N[6]),0.25 *(1+E[6]) * (1 - N[6]),0.25 *(1+E[6]) * (1 + N[6]),0.25 *(1-E[6]) * (1 + N[6])},
                    {0.25 *(1-E[7]) * (1 - N[7]),0.25 *(1+E[7]) * (1 - N[7]),0.25 *(1+E[7]) * (1 + N[7]),0.25 *(1-E[7]) * (1 + N[7])},
                    {0.25 *(1-E[8]) * (1 - N[8]),0.25 *(1+E[8]) * (1 - N[8]),0.25 *(1+E[8]) * (1 + N[8]),0.25 *(1-E[8]) * (1 + N[8])}
            };


        }

    }
    public double[][] solution(int a, int i , List<Element> nElem, List<Node> nNode){
        double ksiX=0;
        double ksiY=0;
        double etaX=0;
        double etaY=0;

        double [] dndx = new double[4];
        double [] dndy = new double[4];
        double [][] resultX= new double[4][4];
        double [][] resultY= new double[4][4];
        double [][] last = new double[4][4];

        for(int j=0;j<4;j++){
            ksiX=ksiX+pktcalkE[i][j]*nNode.get(nElem.get(i).ID[j]-1).x;
            ksiY=ksiY+pktcalkE[i][j]*nNode.get(nElem.get(i).ID[j]-1).y;
            etaX=etaX+pktcalkN[i][j]*nNode.get(nElem.get(i).ID[j]-1).x;
            etaY=etaY+pktcalkN[i][j]*nNode.get(nElem.get(i).ID[j]-1).y;

        }

        double [][] Jackobian = new double [2][2];
        Jackobian[0][0]=ksiX;
        Jackobian[0][1]=ksiY;
        Jackobian[1][0]=etaX;
        Jackobian[1][1]=etaY;

        double detJackobian = Jackobian[0][0]*Jackobian[1][1] - Jackobian[0][1]*Jackobian[1][0];

        double [][] invJackobian = new double [2][2];
        invJackobian[0][0]=(1/detJackobian)*Jackobian[1][1];
        invJackobian[0][1]=(1/detJackobian)*(-Jackobian[0][1]);
        invJackobian[1][0]=(1/detJackobian)*(-Jackobian[1][0]);
        invJackobian[1][1]=(1/detJackobian)*Jackobian[0][0];

        System.out.println("Jackobian dla pkt calkowania nr: " + (i+1));
        for (int k=0;k<2;k++){
            for(int l=0;l<2;l++){
                System.out.print(Jackobian[k][l]+ " ");
            }
            System.out.println();
        }
        System.out.println();
        System.out.println();
        System.out.println("Wyznacznik= " + detJackobian);
        System.out.println();
        System.out.println("Odwrocony Jackobian dla pkt calkowania nr: "+ (i+1));

        for (int n=0;n<2;n++){
            for(int m=0;m<2;m++){
                System.out.print(invJackobian[n][m]+ " ");
            }
            System.out.println();
        }
        System.out.println();

        for (int j=0;j<4;j++){
            dndx[j]=invJackobian[0][0]*pktcalkE[i][j] + invJackobian[0][1]*pktcalkN[i][j];
            dndy[j]=invJackobian[1][0]*pktcalkE[i][j] + invJackobian[1][1]*pktcalkN[i][j];
            System.out.println("dn/dx: " + dndx[j] + " dn/dy: " + dndy[j]);
            System.out.println();
        }

        for (int m=0;m<4;m++){
            for (int n=0; n<4;n++){
                resultX[m][n]=dndx[m]*dndx[n];
                resultY[m][n]=dndy[m]*dndy[n];
                System.out.print(resultX[m][n]);
            }
            System.out.println();
        }
        System.out.println();

        for (int m=0;m<4;m++){
            for (int n=0;n<4;n++){
                System.out.print(resultY[m][n]);
            }
            System.out.println();
        }

        for (int z=0;z<4;z++){
            for (int v=0;v<4;v++){
                last[z][v]=25*detJackobian*(resultX[z][v] + resultY[z][v])*scale1[i]*scale2[i];

            }
        }

        for (int j=0;j<4;j++){
            for(int k=0;k<4;k++){
                System.out.print(last[j][k]+ " ");
            }
            System.out.println();
        }
        System.out.println();

        return last;
    }

    public double[][] solutionC(int a, int i , List<Element> nElem, List<Node> nNode, double p, double cp){
        double  ksiX=0;
        double  ksiY=0;
        double  etaX=0;
        double  etaY=0;
        double [][] var= new double[4][4];
        double [][] C = new double [4][4];

        for(int j=0;j<4;j++){
            ksiX=ksiX+pktcalkE[i][j]*nNode.get(nElem.get(i).ID[j]-1).x;
            ksiY=ksiY+pktcalkE[i][j]*nNode.get(nElem.get(i).ID[j]-1).y;
            etaX=etaX+pktcalkN[i][j]*nNode.get(nElem.get(i).ID[j]-1).x;
            etaY=etaY+pktcalkN[i][j]*nNode.get(nElem.get(i).ID[j]-1).y;

        }

        double [][] Jackobian = new double [2][2];
        Jackobian[0][0]=ksiX;
        Jackobian[0][1]=ksiY;
        Jackobian[1][0]=etaX;
        Jackobian[1][1]=etaY;

        double detJackobian = Jackobian[0][0]*Jackobian[1][1] - Jackobian[0][1]*Jackobian[1][0];

        for(int s=0;s<4;s++)
        {
            for(int d=0;d<4;d++)
            {
                var[s][d]=L[i][s]*L[i][d];

            }
        }

        for(int g=0;g<4;g++){
            for(int h=0;h<4;h++){
                C[g][h] = ((var[g][h])*p*cp*detJackobian)*scale1[i]*scale2[i];
            }
        }
        return C;
    }

    private double calcN (int a, double ksi, double eta ){
        if(a==0)return 0.25*(1-ksi)*(1-eta);
        if(a==1)return 0.25*(1+ksi)*(1-eta);
        if(a==2)return 0.25*(1+ksi)*(1+eta);
        if(a==3)return 0.25*(1-ksi)*(1+eta);
        else {
            return 0;
        }

    }

    public void calcH(List<Element> nElem, List<Node> nodes,double alfa){
        for (Element element : nElem) {
            element.addtoHL(giveHBC(element,nodes,alfa));

        }
    }
    private double [][]giveHBC (Element element,List<Node> nodes,double alfa){

        double [][]HBC  = new double [4][4];
        Node []elemNodes = new Node[]{
            nodes.get(element.ID[0]-1),
            nodes.get(element.ID[1]-1),
            nodes.get(element.ID[2]-1),
            nodes.get(element.ID[3]-1)
        };


        if (numberP==2) {
            double []ksi=new double[2];
            double []eta=new double[2];
            if (elemNodes[0].bc == 1 && elemNodes[1].bc == 1) {
                eta[0] = -1;
                eta[1] = -1;
                ksi[0] = (-1 / Math.sqrt(3));
                ksi[1] = (1 / Math.sqrt(3));
                double l = Math.abs((elemNodes[0].x - elemNodes[1].x) / 2);
                calcSides (ksi, eta, HBC, alfa, l);
            }

            if (elemNodes[1].bc == 1 && elemNodes[2].bc == 1) {
                eta[0] = (-1 / Math.sqrt(3));
                eta[1] = (1 / Math.sqrt(3));
                ksi[0] = 1;
                ksi[1] = 1;
                double l = Math.abs((elemNodes[1].y - elemNodes[2].y) / 2);
                calcSides(ksi, eta, HBC, alfa, l);
            }

            if (elemNodes[2].bc == 1 && elemNodes[3].bc == 1) {
                eta[0] = 1;
                eta[1] = 1;
                ksi[0] = (-1 / Math.sqrt(3));
                ksi[1] = (1 / Math.sqrt(3));
                double l = Math.abs((elemNodes[2].x - elemNodes[3].x) / 2);
                calcSides(ksi, eta, HBC, alfa, l);
            }

            if (elemNodes[3].bc == 1 && elemNodes[0].bc == 1) {
                eta[0] = (-1 / Math.sqrt(3));
                eta[1] = (1 / Math.sqrt(3));
                ksi[0] = -1;
                ksi[1] = -1;
                double l = Math.abs((elemNodes[3].y - elemNodes[0].y) / 2);
                calcSides(ksi, eta, HBC, alfa, l);
            }

        }
        else if (numberP==3) {
            double []ksi=new double[3];
            double []eta=new double[3];
            if (elemNodes[0].bc == 1 && elemNodes[1].bc == 1) {
                eta[0] = -1;
                eta[1] = -1;
                eta[2] = -1;
                ksi[0] = (-1 * Math.sqrt(3.0/5.0));
                ksi[1] = 0;
                ksi[2] = (1 * Math.sqrt(3.0/5.0));
                double l = Math.abs((elemNodes[0].x - elemNodes[1].x) / 2);
                calcSides(ksi, eta, HBC, alfa, l);
            }

            if (elemNodes[1].bc == 1 && elemNodes[2].bc == 1) {
                eta[0] = (-1 * Math.sqrt(3.0/5.0));
                eta[1] = 0;
                eta[2] = (1 * Math.sqrt(3.0/5.0));
                ksi[0] = 1;
                ksi[1] = 1;
                ksi[2] = 1;
                double l = Math.abs((elemNodes[1].y - elemNodes[2].y) / 2);
                calcSides(ksi, eta, HBC, alfa, l);
            }

            if (elemNodes[2].bc == 1 && elemNodes[3].bc == 1) {
                eta[0] = 1;
                eta[1] = 1;
                eta[2] =1;
                ksi[0] = (-1 * Math.sqrt(3.0/5.0));
                ksi[1] = 0;
                ksi[2] =(1 * Math.sqrt(3.0/5.0));
                double l = Math.abs((elemNodes[2].x - elemNodes[3].x) / 2);
                calcSides(ksi, eta, HBC, alfa, l);
            }

            if (elemNodes[3].bc == 1 && elemNodes[0].bc == 1) {
                eta[0] = (-1 * Math.sqrt(3.0/5.0));
                eta[1] = 0;
                eta[2]=(1 * Math.sqrt(3.0/5.0));
                ksi[0] = -1;
                ksi[1] = -1;
                ksi[2]=-1;
                double l = Math.abs((elemNodes[3].y - elemNodes[0].y) / 2);
                calcSides(ksi, eta, HBC, alfa, l);
            }

        }
        return HBC;

    }
    private void calcSides (double []ksi, double []eta, double [][]HBC, double alfa, double length ){

        double [][][] sidePointsH = new double [numberP][4][4];
        double [][] sideVector = new double [numberP][4];
        double [][] sideHBC = new double [4][4];
        for (int i = 0; i < numberP; i++) {
            for (int i1 = 0; i1 < 4; i1++) {
                sideVector[i][i1]=calcN(i1,ksi[i],eta[i]);

            }
        }

        for (int i = 0; i < numberP; i++) {
            for (int j = 0; j <4 ; j++) {
                for (int k = 0; k < 4; k++) {
                    sidePointsH[i][j][k]=sideVector[i][j]*sideVector[i][k];
                }
            }
        }

        for (int i = 0; i < numberP; i++) {
            for (int j = 0; j <4 ; j++) {
                for (int k = 0; k < 4; k++) {
                    sideHBC[j][k]+=sidePointsH[i][j][k]*scale1[i];
                }
            }
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                HBC[i][j]+=sideHBC[i][j]*alfa*length;
            }
        }

    }

    public void calcPC(List<Element> nElem, List<Node> nodes,double alfa, double area){
        for (Element element : nElem) {
            element.setPC(givePC(element,nodes,alfa,area));

        }
    }
    private double []givePC (Element element,List<Node> nodes,double alfa, double area) {
        double[] PC = new double[4];
        Node[] elemNodes = new Node[]{
                nodes.get(element.ID[0] - 1),
                nodes.get(element.ID[1] - 1),
                nodes.get(element.ID[2] - 1),
                nodes.get(element.ID[3] - 1)
        };
        if (numberP == 2) {
            double[] ksi = new double[2];
            double[] eta = new double[2];

            if (elemNodes[0].bc == 1 && elemNodes[1].bc == 1) {
                eta[0] = -1;
                eta[1] = -1;
                ksi[0] = (-1 / Math.sqrt(3));
                ksi[1] = (1 / Math.sqrt(3));
                double l = Math.abs((elemNodes[0].x - elemNodes[1].x) / 2);
                calcSidesPC(ksi, eta, PC, alfa, l, area);
            }

            if (elemNodes[1].bc == 1 && elemNodes[2].bc == 1) {
                eta[0] = (-1 / Math.sqrt(3));
                eta[1] = (1 / Math.sqrt(3));
                ksi[0] = 1;
                ksi[1] = 1;
                double l = Math.abs((elemNodes[1].y - elemNodes[2].y) / 2);
                calcSidesPC(ksi, eta, PC, alfa, l, area);
            }

            if (elemNodes[2].bc == 1 && elemNodes[3].bc == 1) {
                eta[0] = 1;
                eta[1] = 1;
                ksi[0] = (-1 / Math.sqrt(3));
                ksi[1] = (1 / Math.sqrt(3));
                double l = Math.abs((elemNodes[2].x - elemNodes[3].x) / 2);
                calcSidesPC(ksi, eta, PC, alfa, l, area);
            }

            if (elemNodes[3].bc == 1 && elemNodes[0].bc == 1) {
                eta[0] = (-1 / Math.sqrt(3));
                eta[1] = (1 / Math.sqrt(3));
                ksi[0] = -1;
                ksi[1] = -1;
                double l = Math.abs((elemNodes[3].y - elemNodes[0].y) / 2);
                calcSidesPC(ksi, eta, PC, alfa, l, area);
            }
        } else if (numberP == 3) {
            double[] ksi = new double[3];
            double[] eta = new double[3];

            if (elemNodes[0].bc == 1 && elemNodes[1].bc == 1) {
                eta[0] = -1;
                eta[1] = -1;
                eta[2] = -1;
                ksi[0] = (-1 * Math.sqrt(3.0/5.0));
                ksi[1] = 0;
                ksi[2] = (1 * Math.sqrt(3.0/5.0));
                double l = Math.abs((elemNodes[0].x - elemNodes[1].x) / 2);
                calcSidesPC(ksi, eta, PC, alfa, l, area);
            }

            if (elemNodes[1].bc == 1 && elemNodes[2].bc == 1) {
                eta[0] = (-1 * Math.sqrt(3.0/5.0));
                eta[1] = 0;
                eta[2] = (1 * Math.sqrt(3.0/5.0));
                ksi[0] = 1;
                ksi[1] = 1;
                ksi[2] = 1;
                double l = Math.abs((elemNodes[1].y - elemNodes[2].y) / 2);
                calcSidesPC(ksi, eta, PC, alfa, l, area);
            }

            if (elemNodes[2].bc == 1 && elemNodes[3].bc == 1) {
                eta[0] = 1;
                eta[1] = 1;
                eta[2] =1;
                ksi[0] = (-1 * Math.sqrt(3.0/5.0));
                ksi[1] = 0;
                ksi[2] =(1 * Math.sqrt(3.0/5.0));
                double l = Math.abs((elemNodes[2].x - elemNodes[3].x) / 2);
                calcSidesPC(ksi, eta, PC, alfa, l, area);
            }

            if (elemNodes[3].bc == 1 && elemNodes[0].bc == 1) {
                eta[0] = (-1 * Math.sqrt(3.0/5.0));
                eta[1] = 0;
                eta[2] = (1 * Math.sqrt(3.0/5.0));
                ksi[0] = -1;
                ksi[1] = -1;
                ksi[2] =-1;
                double l = Math.abs((elemNodes[3].y - elemNodes[0].y) / 2);
                calcSidesPC(ksi, eta, PC, alfa, l, area);

            }

        }
        return PC;
    }
    private void calcSidesPC (double []ksi, double []eta, double []PC, double alfa, double l,double area ){
        double [] sidePC = new double [4];
        double [][] pointPC = new double [4][4];
        for (int i = 0; i < numberP; i++) {
            for (int i1 = 0; i1 < 4; i1++) {
                pointPC[i][i1]=calcN(i1,ksi[i],eta[i]);

            }
        }

        for (int i = 0; i < numberP; i++) {
            for (int j = 0; j <4 ; j++) {
                    sidePC[j]+=pointPC[i][j]*scale1[i];
            }
        }


        for (int j = 0; j < 4; j++) {
            PC[j]-=sidePC[j]*alfa*area*l;
        }


    }



}
