import java.util.ArrayList;
import java.util.List;


public class Main {

    public static void main(String[] args) {
        Global global = new Global();
        int pc= global.pc;
        double dX=global.W/(global.nW-1);
        double dY=global.H/(global.nH-1);
        double nE=global.nE;
        double nN=global.nN;
        double [][] temps;
        double [][] tempsC;
       List<Element>nElem=new ArrayList<>();
       List<Node>nNode=new ArrayList<>();

       for(int i=0; i<global.nW; i++){
           for (int j=0; j<global.nH;j++){
               double x =i*dX;
               double y =j*dY;
               nNode.add(new Node(x,y));
           }
       }
       for(int i=1;i<global.nE+(global.nW-1); i++){
           if (i%global.nH==0){
               continue;
           }
           Element tempElem = new Element();
           tempElem.ID[0]=i;
           tempElem.ID[1]=tempElem.ID[0]+global.nH;
           tempElem.ID[2]=tempElem.ID[1]+1;
           tempElem.ID[3]=tempElem.ID[0]+1;
           nElem.add(tempElem);
       }

       IElem3 matrix;
       if(pc==2){
            matrix = new Element3(2);
       }
       else if(pc==3) {
           matrix = new Element3(3);
       }
       else{
           matrix = new Element3(2);
       }
       double [][] HG = new double [global.nN][global.nN];
       matrix.calcH(nElem,nNode, global.alfa);
       for( int p=0;p<global.nE;p++){
           for(int i=0;i< global.pc* global.pc;i++){
               temps=matrix.solution(p,i,nElem,nNode);
               for(int k=0;k<4;k++){
                   for (int l=0;l<4;l++){
                        nElem.get(p).Hl[k][l] += temps[k][l];
                   }
               }
           }
           System.out.println();
           System.out.println("Hl Matrix");
           for(int k=0;k<4;k++){
               for(int j=0;j<4;j++){
                   System.out.print(nElem.get(p).Hl[k][j]+ " ");
               }
               System.out.println();
           }
       }
           System.out.println();
           for(int t=0;t<global.nE;t++){
               for (int i=0;i<4;i++){
                   for(int j=0;j<4;j++){
                       HG[nElem.get(t).ID[i]-1][nElem.get(t).ID[j]-1]+=nElem.get(t).Hl[i][j];
                   }
               }

           }

           for (int i=0;i< global.nN; i++){
               for (int k=0;k< global.nN;k++){
                   System.out.print(HG[i][k]+" ");
               }
               System.out.println();
           }


        double [][] CG=new double [global.nN][global.nN];
        for( int p=0;p<global.nE;p++){
            for(int i=0;i<global.pc* global.pc;i++){
                tempsC=matrix.solutionC(p,i,nElem,nNode,global.p, global.cp);
                for(int k=0;k<4;k++){
                    for (int l=0;l<4;l++){
                        nElem.get(p).Cl[k][l] += tempsC[k][l];
                    }
                }
            }
            System.out.println();
            System.out.println("Cl Matrix");
            for(int k=0;k<4;k++){
                for(int j=0;j<4;j++){
                    System.out.print(nElem.get(p).Cl[k][j]+ " ");

                }
                System.out.println();
            }
        }
        System.out.println();
        for(int t=0;t<global.nE;t++){
            for (int i=0;i<4;i++){
                for(int j=0;j<4;j++){
                    CG[nElem.get(t).ID[i]-1][nElem.get(t).ID[j]-1]+=nElem.get(t).Cl[i][j];
                }
            }

        }

        for (int i=0;i< global.nN; i++){
            for (int k=0;k< global.nN;k++){
                System.out.print(CG[i][k]+" ");
            }
            System.out.println();
        }

        matrix.calcPC(nElem,nNode, global.alfa, global.otoczenie);
        for (Element element : nElem) {
            for (int i = 0; i < element.PC.length; i++) {
                System.out.print(element.PC[i]+" ");
            }
            System.out.println();
        }

        double []PG = new double [nNode.size()];

        for (int i = 0; i < nElem.size() ; i++) {
            for (int j = 0; j < 4; j++) {
                PG[(nElem.get(i).ID[j]-1)]+=nElem.get(i).PC[j];
            }
        }
        Simulation sim = new Simulation();
        sim.count(HG,CG,PG,global);

       }



}

