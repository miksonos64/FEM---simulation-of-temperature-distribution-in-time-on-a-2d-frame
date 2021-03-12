import java.util.List;

public interface IElem3 {
    public double[][] solution(int a, int i , List<Element> nElem, List<Node> nNode);
    public double[][] solutionC(int a, int i , List<Element> nElem, List<Node> nNode, double p, double cp);
    public void calcH(List<Element> nElem, List<Node> nodes,double alfa);
    public void calcPC(List<Element> nElem, List<Node> nodes,double alfa, double area);
}
