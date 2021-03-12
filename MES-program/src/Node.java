public class Node {
    public double x, y;
    public double t0;
    public int bc;
    public Node(double x, double y){
        this.x = x;
        this.y = y;
        this.t0=t0;
        if(this.x ==0 || this.y==0 || this.x == new Global().H || this.y == new Global().W ){
            this.bc =1;
        }
        else{
            this.bc=0;
        }
    }
}