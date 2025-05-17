package adduct;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/* Adduct: product of a direct addition of two or more distinct molecules,
 * resulting in a single reaction product containing all atoms of all components.
 */
public abstract class Adduct {

    /**
     * Calculate the mass to search depending on the adduct hypothesis.
     * Supports single and multiple charges, as well as multimers (e.g., [2M+Na]+).
     *
     * @param mz mz
     * @param adduct adduct notation ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, Map.Entry<String, Double> adduct) {

        Double monoisotopicMass = null;
        // !! TODO METHOD
        // !! TODO Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).
        int multimer = extractMultimer(adduct.getKey());
        int charge = extractCharge(adduct.getKey());

        // If charge is higher than one the adduct shift has to be divided by the charge
        Double adductShift = adduct.getValue()/charge;

        // Case 1: Single charge, no multimer (m/z = M +- adductShift).
        if (charge == 1 && multimer == 1) {
            monoisotopicMass = mz + adductShift;
        }
        // Case 2: Multiple charges, no multimer (mz = M/charge +- adductShift).
        else if (charge > 1 && multimer == 1) {
            monoisotopicMass = (mz + adductShift) * charge;
        }
        // Case 3: Multimer, single charge (mz = M * numberOfMultimer +- adductShift).
        else if (charge == 1 && multimer > 1) {
            monoisotopicMass = (mz + adductShift) / multimer;
        }
        // Case 4: Multimer with multiple charges (mz = (M * numberOfMultimer +- adductShift) / charge).
        else {
            monoisotopicMass = ((mz + adductShift) * charge) / multimer;
        }
        return monoisotopicMass;
    }

    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     * Supports single and multiple charges, as well as multimers (e.g., [2M+Na]+).
     *
     * @param monoisotopicMass The neutral mass of the molecule (M).
     * @param adduct adduct notation (e.g. "[M+H]+", "[2M+Na]+", "[M+2H]2+").
     *
     * @return The expected m/z for the adduct.
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, Map.Entry<String, Double> adduct) {

        Double mz = null;
        // !! TODO METHOD
        // !! TODO Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).
        int multimer = extractMultimer(adduct.getKey());
        int charge = extractCharge(adduct.getKey());

        // If charge is higher than one the adduct shift has to be divided by the charge
        Double adductShift = adduct.getValue()/charge;

        // Case 1: Single charge, no multimer (m/z = M +- adductShift).
        if (charge == 1 && multimer == 1) {
            mz = monoisotopicMass - adductShift;
        }
        // Case 2: Multiple charges, no multimer (mz = M/charge +- adductShift).
        else if (charge > 1 && multimer == 1) {
            mz = (monoisotopicMass / charge) - adductShift;
        }
        // Case 3: Multimer, single charge (mz = M * numberOfMultimer +- adductShift).
        else if (charge == 1 && multimer > 1) {
            mz = (monoisotopicMass * multimer) - adductShift;
        }
        // Case 4: Multimer with multiple charges (mz = (M * numberOfMultimer +- adductShift) / charge).
        else {
            mz = ((monoisotopicMass * multimer)/charge) - adductShift;
        }
        return mz;
    }

    /**
     * Calculates the parts-per-million (PPM) difference between an experimental mass and a theoretical mass.
     *
     * @param experimentalMass the measured or observed mass value
     * @param theoreticalMass the expected or reference mass value
     * @return the absolute PPM difference between the two mass values, rounded to the nearest integer
     */
    public static int calculatePPMIncrement(Double experimentalMass, Double theoreticalMass) {
        int ppmIncrement;
        ppmIncrement = (int) Math.round(Math.abs((experimentalMass - theoreticalMass)*1000000 / theoreticalMass));
        return ppmIncrement;
    }

    /**
     * Calculates the absolute delta in mass corresponding to a given parts-per-million (PPM) tolerance.
     *
     * @param experimentalMass the experimental mass (in Daltons or any relevant mass unit)
     * @param ppm the parts-per-million tolerance
     * @return the absolute mass difference (delta) corresponding to the given PPM, rounded to the nearest integer
     */
    public static double calculateDeltaPPM(Double experimentalMass, int ppm) {
        double deltaPPM;
        deltaPPM = Math.round(Math.abs((experimentalMass * ppm) / 1000000));
        return deltaPPM;
    }

    /**
     * Extracts the multimer count from an adduct string in the format commonly used in mass spectrometry.
     *
     * @param adduct the adduct string
     * @return the integer multiplier before M
     */
    private static int extractMultimer(String adduct) {

        // Pattern looks for a number before 'M', e.g. [2M or [M
        Matcher m = Pattern.compile("\\[([0-9]*)M").matcher(adduct);
        if (m.find()) {
            String num = m.group(1); // Group 1 captures the digits before 'M'
            return num.isEmpty() ? 1 : Integer.parseInt(num);
        }
        return 1; // Default: no multimer specified
    }

    /**
     * Extracts the ionic charge from an adduct string typically used in mass spectrometry.
     * The charge is found as a number directly before the '+' or '-' at the end of the adduct.
     * If no number is present, the method assumes a default charge of 1.
     *
     * @param adduct The adduct string in standard format
     * @return The charge value
     */
    private static int extractCharge(String adduct) {

        // Match the last digit(s) followed by + or − before the closing bracket
        Matcher m = Pattern.compile("([0-9]*)([+-])\\]?$").matcher(adduct);
        if (m.find()) {
            String num = m.group(1);  // May be empty
            return num.isEmpty() ? 1 : Integer.parseInt(num);
        }
        return 1; // Default if no explicit charge found
    }
}
