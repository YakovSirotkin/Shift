package ru.spbau.bioinf.shift;

public class ProteinPosition {
    private int proteinId;
    private double pos;

    public ProteinPosition(int proteinId, double pos) {
        this.proteinId = proteinId;
        this.pos = pos;
    }

    public int getProteinId() {
        return proteinId;
    }

    public double getPos() {
        return pos;
    }
}
