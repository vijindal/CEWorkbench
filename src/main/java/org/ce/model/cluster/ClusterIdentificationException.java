package org.ce.model.cluster;

/**
 * Thrown when a precondition of a {@link ClusterCFIdentificationPipeline} stage
 * is not met (missing/malformed inputs, degenerate intermediate results).
 *
 * <p>Distinguishes expected, diagnosable input problems from programming bugs
 * (NPE, ArrayIndexOutOfBoundsException, etc.) so callers can surface a clear
 * message to the user instead of a raw stack trace.</p>
 */
public class ClusterIdentificationException extends RuntimeException {

    public ClusterIdentificationException(String stage, String message) {
        super("[" + stage + "] " + message);
    }

    public ClusterIdentificationException(String stage, String message, Throwable cause) {
        super("[" + stage + "] " + message, cause);
    }
}
