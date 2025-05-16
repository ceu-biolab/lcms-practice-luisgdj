package lipid;

import adduct.Adduct;
import adduct.AdductList;
import java.util.*;

/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private final IoniationMode ionizationMode;
    private String adduct; // !!TODO The adduct will be detected based on the groupedSignals
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode) {
        this(lipid, mz, intensity, retentionTime, ionizationMode, Collections.emptySet());
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionizationMode = ionizationMode;
        // !!TODO This set should be sorted according to help the program to deisotope the signals plus detect the adduct
        this.groupedSignals = new TreeSet<>(Comparator.comparing(Peak::getMz)); // Sorted based on m/z value
        this.groupedSignals.addAll(groupedSignals);
        this.score = 0;
        this.totalScoresApplied = 0;
        detectAdductFromPeaks();
    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public IoniationMode getIonizationMode() {
        return ionizationMode;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }

    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    // !TODO Take into account that the score should be normalized between 0 and 1
    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    public double getNormalizedScore() {
        if (this.totalScoresApplied == 0) {
            return 0; // or Double.NaN or throw exception, depending on desired behavior
        } else {
            return (double) this.score / this.totalScoresApplied;
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 && Double.compare(that.rtMin, rtMin) == 0 && Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }

    // !!TODO Detect the adduct with an algorithm or with drools, up to the user.
    /**
     * Tries to detect the adduct type for the current annotation based on grouped MS peaks.
     * It compares all pairs of peaks, assuming each peak is an adduct form of the same neutral molecule,
     * and checks if their theoretical mass difference fits known adduct transformations.
     */
    private void detectAdductFromPeaks() {

        List<Peak> peaks = new ArrayList<>(groupedSignals);

        for (int i = 0; i < peaks.size(); i++) {
            Peak p1 = peaks.get(i);

            if(Math.abs(p1.getMz() - this.mz) < 0.01) {
                for (int j = 0; j < peaks.size(); j++) {
                    Peak p2 = peaks.get(j);
                    if(p1.equals(p2)) continue;
                    System.out.println("Peak 1: " + p1 + "\nPeak 2: " + p2);

                    if(ionizationMode.equals(IoniationMode.POSITIVE)) {
                        // Try detecting a positive adduct match
                        this.adduct = detectAdductFromMz(p1, p2, AdductList.MAPMZPOSITIVEADDUCTS);
                        return;
                    }
                    if(ionizationMode.equals(IoniationMode.NEGATIVE)) {
                        // Try detecting a positive adduct match
                        this.adduct = detectAdductFromMz(p1, p2, AdductList.MAPMZNEGATIVEADDUCTS);
                        return;
                    }
                }
            }
        }
        System.out.println("ERROR: No adduct detected for annotation with mz = " + this.mz);
    }

    public static String detectAdductFromMz(Peak p1, Peak p2, Map<String,Double> mapMz) {

        // Try every combination of known positive ion adducts for both peaks
        for (Map.Entry<String,Double> adduct1 : mapMz.entrySet()) {
            String a1 = adduct1.getKey();
            Double mz1 = p1.getMz();
            // Predict the monoisotopic mass value of the first peak with its m/z and adduct a1
            Double theoreticalM1 = Adduct.getMonoisotopicMassFromMZ(mz1, adduct1);

            for (Map.Entry<String,Double> adduct2 : mapMz.entrySet()) {
                String a2 = adduct2.getKey();
                Double mz2 = p2.getMz();
                // Predict the m/z value of the second peak assuming the same neutral mass M1 with adduct a2
                Double theoreticalMz2 = Adduct.getMZFromMonoisotopicMass(theoreticalM1, adduct2);

                // If the adduct types are different, evaluate how close the measured and theoretical m/z are
                if (!a1.equals(a2)) {
                    // Define the accepted error tolerance in parts per million (ppm)
                    int ppmTolerance = 10;
                    int ppmError = Adduct.calculatePPMIncrement(mz2, theoreticalMz2);
                    // If the ppm error is within tolerance, it's a match
                    if (ppmError <= ppmTolerance) {
                        // Set the adduct type of the annotation based on the first peak
                        System.out.println("Detected pair:");
                        System.out.println(" - Adducts: " + a1 + " & " + a2);
                        System.out.println(" - Theoretical M: " + theoreticalM1);
                        return a1;
                    }
                }
            }
        }
        return null;
    }
}
