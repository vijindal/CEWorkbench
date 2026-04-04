package org.ce.domain.cluster;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Mapping from site-operator products to CF indices.
 */
public final class SubstituteRules {

    private final Map<SiteOpProductKey, CFIndex> rules;

    SubstituteRules(Map<SiteOpProductKey, CFIndex> rules) {
        this.rules = Collections.unmodifiableMap(new LinkedHashMap<>(rules));
    }

    public CFIndex lookup(List<SiteOp> ops) {
        return rules.get(new SiteOpProductKey(ops));
    }

    public Map<SiteOpProductKey, CFIndex> getRules() {
        return rules;
    }
}

