package org.ce.storage;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.function.Function;

/**
 * Generic storage class for JSON-serialized data.
 * 
 * <p>Parameterized by the type of data being stored and a path resolver 
 * function that maps an ID to a filesystem {@link Path}.</p>
 */
public class DataStore<T> {

    private final Function<String, Path> pathResolver;
    private final Class<T> type;
    private final ObjectMapper mapper = new ObjectMapper();

    public DataStore(Function<String, Path> pathResolver, Class<T> type) {
        this.pathResolver = pathResolver;
        this.type = type;
    }

    public boolean exists(String id) {
        return Files.exists(pathResolver.apply(id));
    }

    public T load(String id) throws IOException {
        Path file = pathResolver.apply(id);
        if (!Files.exists(file)) {
            throw new IOException(type.getSimpleName() + " not found: " + file);
        }
        return mapper.readValue(file.toFile(), type);
    }

    public void save(String id, T data) throws IOException {
        Path file = pathResolver.apply(id);
        Files.createDirectories(file.getParent());
        mapper.writerWithDefaultPrettyPrinter()
              .writeValue(file.toFile(), data);
    }
}
