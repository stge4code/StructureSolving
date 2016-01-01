package CrystTools;

/**
 * Created by Developer on 10.09.2015.
 */
public class ReciprocalItem {
    public int h;
    public int k;
    public int l;
    public double Fsq;
    public double sigmaFsq;
    public double batchNumber;
    public double scatvect;

    public ReciprocalItem(int h, int k, int l, double fsq, double sigmaFsq, double scatvect, double batchnumber) {
        this.h = h;
        this.k = k;
        this.l = l;
        this.Fsq = fsq;
        this.sigmaFsq = sigmaFsq;
        this.scatvect = scatvect;
        this.batchNumber = batchnumber;
    }
    public ReciprocalItem(ReciprocalItem other) {
        this.h = other.h;
        this.k = other.k;
        this.l = other.l;
        this.Fsq = other.Fsq;
        this.sigmaFsq = other.sigmaFsq;
        this.scatvect = other.scatvect;
        this.batchNumber = other.batchNumber;
    }
}
