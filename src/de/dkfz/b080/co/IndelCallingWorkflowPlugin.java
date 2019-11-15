/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
 *
 */
package de.dkfz.b080.co;

import de.dkfz.roddy.plugins.BasePlugin;

/**
 */
public class IndelCallingWorkflowPlugin extends BasePlugin {

    public static final String CURRENT_VERSION_STRING = "2.4.1";
    public static final String CURRENT_VERSION_BUILD_DATE = "Fri Nov 15 12:41:48 CET 2019";

    @Override
    public String getVersionInfo() {
        return "Roddy plugin: " + this.getClass().getName() + ", V " + CURRENT_VERSION_STRING + " built at " + CURRENT_VERSION_BUILD_DATE;
    }
}
