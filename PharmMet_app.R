library(shiny)
library(shinythemes)
library(DT)
library(tidyr)
library(dplyr)
library(fontawesome)
library(shinyBS)

# Read in dataset
data <- read.csv("PharmMet_DB_061325.csv")

normalize_adduct <- function(x) {
  x <- gsub("'", "", x)             
  x <- gsub("\\[|\\]", "", x)        
  x <- trimws(x)                     
  x <- sub("([+-])$", "", x)         
  return(x)
}

adduct_choices_raw <- c("[M+2H]+", "[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+",
                        "[2M+H]+", "[M-2H]2-", "[M-H2O-H]-", "[M-H]-", "[M+Na-2H]-", "[M+Cl]-")

adduct_choices <- normalize_adduct(adduct_choices_raw)

expand_dataset <- function(data) {
  data <- data %>%
    mutate(
      adduct = normalize_adduct(adduct),
      mz = gsub("\\[|\\]", "", mz),
      rt = gsub("\\[|\\]", "", rt),
      mode = gsub("'", "", gsub("\\[|\\]", "", mode))
    ) %>%
    separate_rows(adduct, mz, rt, mode, sep = ", ") %>%
    mutate(
      adduct = trimws(adduct),
      mz = suppressWarnings(as.numeric(trimws(mz))),
      rt = suppressWarnings(as.numeric(trimws(rt)))  
    ) %>%
    mutate(
      rt = round(rt * 60, 0)  
    )
  
  if("Reaction_phases" %in% names(data)) {
    data <- data %>% 
      mutate(Reaction_phases = ifelse(is.na(enzymes) | enzymes == "", "Unknown", Reaction_phases))
  } else {
    data <- data %>% mutate(Reaction_phases = ifelse(is.na(enzymes) | enzymes == "", "Unknown", "Available"))
  }
  
  data <- data %>% 
    mutate(detected = ifelse(is.na(detected) | detected == "", "Not detected", detected)) %>%
    filter(!is.na(mz), !is.na(rt)) %>%
    distinct(name, mz, rt, adduct, mode, detected, .keep_all = TRUE)
  
  return(data)
}

# Expanded dataset
expanded_data <- expand_dataset(data)

# UI
ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com"),
    tags$link(
      rel = "stylesheet", 
      href = "https://fonts.googleapis.com/css2?family=Montserrat:wght@400;500;700&display=swap"
    ),
    tags$style(shiny::HTML("
      body {
        font-family: 'Montserrat', sans-serif;
        background-color: #f9fbfc;
        color: #343a40;
      }
      .navbar { background: linear-gradient(to right, #2c3e50, #34495e) !important; border: none; }
      .navbar-brand { font-weight: 700; color: #ffffff !important; }
      .navbar-nav > li > a { color: #ffffff !important; font-size: 16px; font-weight: 500; }
      .main-header { background: linear-gradient(to right, #2980b9, #3498db); padding: 30px;
                     text-align: center; border-radius: 8px; margin-bottom: 30px; }
      .main-header h1 { color: white !important; font-size: 48px; margin-bottom: 10px; font-weight: 700; }
      .main-header p { color: white !important; font-size: 20px; font-weight: 400; }
      h3 { font-family: 'Montserrat', sans-serif; font-weight: 700; color: #2c3e50; }
      .card { background-color: #ffffff; border-radius: 10px; padding: 25px;
              box-shadow: 0px 2px 10px rgba(0,0,0,0.05); margin-bottom: 20px;
              transition: box-shadow 0.3s ease; }
      .card:hover { box-shadow: 0px 5px 20px rgba(0,0,0,0.1); }
      .btn-custom { background-color: #1abc9c !important; color: #ffffff !important;
                    border: none; border-radius: 8px; font-size: 16px; padding: 10px 20px;
                    margin-top: 15px; width: 100%; font-weight: 500; transition: background-color 0.3s ease; }
      .btn-custom:hover { background-color: #16a085 !important; }
      .dataTables_wrapper .dataTables_paginate .paginate_button { background-color: #1abc9c;
            color: white !important; border-radius: 4px; margin: 2px; }
      .dataTables_wrapper .dataTables_paginate .paginate_button:hover { background-color: #16a085; }
      tr.details-shown { background-color: #f7f9fa !important; }
      tr.details-shown td, tr.details-shown div { background-color: #f7f9fa !important;
                                                  padding: 10px; border-left: 4px solid #1abc9c; }
      table.dataTable tbody tr:hover { background-color: #f1f5f9 !important; cursor: pointer;
                                       transition: background-color 0.3s ease; }
      footer { padding: 20px; background: linear-gradient(to right, #2c3e50, #34495e);
              color: white; text-align: center; margin-top: 50px; border-radius: 8px; }
      .fa.flask { margin-right: 8px; }
      .header-logo { margin-right: 10px; }
    "))
  ),
  navbarPage(
    title = div(icon("flask", class = "fa-lg header-logo"), "PharmMet Database"),
    
    # Drug-Based Search Tab
    tabPanel("Drug-Based Search",
             div(class = "main-header",
                 h1(icon("flask"), "PharmMet Database"),
                 p("*Pharmaco-metabolomics database for high throughput drug detection*")
             ),
             fluidRow(
               column(3,
                      div(class = "card",
                          h3("Drug-Based Search"),
                          checkboxInput("search_drug_name", "Search by Drug Name", value = TRUE),
                          conditionalPanel(
                            condition = "input.search_drug_name == true",
                            textInput("drug_name", label = icon("search"), placeholder = "Enter Drug Name")
                          ),
                          bsCollapse(id = "drug_group_collapse", open = FALSE, multiple = FALSE,
                                     bsCollapsePanel("Drug Group Selection",
                                                     checkboxGroupInput("drug_group", "Select Drug Group",
                                                                        choices = c("Antiinflammatory and Immunomodulatory", "Antimicrobial",
                                                                                    "Cardiovascular", "Central nervous system", "Endocrine and metabolic",
                                                                                    "Gastrointestinal", "Hematological",
                                                                                    "Musculoskeletal", "Respiratory", "Others"),
                                                                        selected = NULL)
                                     )
                          ),
                          bsCollapse(id = "drug_class_collapse", open = FALSE, multiple = FALSE,
                                     bsCollapsePanel("Drug Class Selection",
                                                     checkboxGroupInput("drug_class", "Select Drug Class",
                                                                        choices = c("Anti-allergic", "Anti-diabetic", "Anti-hypertensive",
                                                                                    "Anti-inflammatory", "Anti-obesity", "Anti-psychotic", "Anti-viral"),
                                                                        selected = NULL)
                                     )
                          ),
                          bsCollapse(id = "reaction_phase_collapse", open = FALSE, multiple = FALSE,
                                     bsCollapsePanel("Reaction Phase Selection",
                                                     checkboxGroupInput("phase_selection", "Select Reaction Phase",
                                                                        choices = c("Phase I", "Phase II", "Unknown"),
                                                                        selected = NULL)
                                     )
                          ),
                          bsCollapse(id = "reaction_type_collapse", open = FALSE, multiple = FALSE,
                                     bsCollapsePanel("Reaction Type Selection",
                                                     checkboxGroupInput("reaction_simplified", "Select Reaction Type",
                                                                        choices = c("Oxidation", "Reduction", "Hydrolysis", "Acetylation",
                                                                                    "Glucuronidation", "Amino acid conjugation", "Carnitine conjugation", "Sulfation",
                                                                                    "GSH-conjugation", "Methylation", "Carboxylation", "Glycosylation","Unknown"),
                                                                        selected = NULL)
                                     )
                          ),
                          actionButton("drug_search_button", label = "Search", class = "btn btn-custom")
                      )
               ),
               column(9,
                      div(class = "card",
                          h3("Search Results"),
                          DTOutput("drug_table")
                      )
               )
             )
    ),
    
    # LC-MS Based Search Tab
    tabPanel("LC-MS Based Search",
             div(class = "main-header",
                 h1(icon("flask"), "PharmMet Database"),
                 p("*Pharmaco-metabolomics database for high throughput drug detection*")
             ),
             fluidRow(
               column(3,
                      div(class = "card",
                          h3("LC-MS Based Search"),
                          numericInput("mz_value", "m/z (3 d.p)", value = NA, min = 0, step = 0.0001),
                          numericInput("rt_value", "Retention Time (sec)", value = NA, min = 0, step = 1),
                          numericInput("ppm_tolerance", "m/z Tolerance (ppm)", value = 10, min = 0, step = 1),
                          checkboxGroupInput("adduct_search", "Select Adduct",
                                             choices = adduct_choices,
                                             selected = NULL),
                          actionButton("lcms_search_button", label = "Search", class = "btn btn-custom")
                      )
               ),
               column(9,
                      div(class = "card",
                          h3("Search Results"),
                          DTOutput("lcms_table")
                      )
               )
             )
    ),
    
    # About Tab
    tabPanel("About",
             div(class = "main-header",
                 h3(icon("info-circle"), "About PharmMet Database")
             ),
             div(class = "card",
                 p("This database provides an interactive way to explore metabolomics data related to drugs. 
                   Human liver S9 system was utilized to generate drug metabolites.
                   References for Biotransformer 3.0 and Drugbank information are available."),
                 br(),
                 p("For citations:"),
                 tags$ul(
                   tags$li(strong("Dean P. Jones"), " / ", strong("Young-Mi Go Lab")),
                   tags$li("PharmMet Database"),
                   tags$li("Citation: xxxx")
                 )
             )
    )
  ),
  tags$footer(
    p("© 2024 PharmMet Database | Developed by Dean P. Jones / Young-Mi Go Lab")
  )
)

# Define the server function
server <- function(input, output, session) {
  
  # Drug-Based Search: Filter Data
  filtered_data_drug <- eventReactive(input$drug_search_button, {
    filtered <- expanded_data
    if (input$search_drug_name && nzchar(input$drug_name)) {
      filtered <- filtered %>% filter(grepl(input$drug_name, name, ignore.case = TRUE))
    }
    if (!is.null(input$drug_group) && length(input$drug_group) > 0) {
      filtered <- filtered %>% filter(drug_group %in% input$drug_group)
    }
    if (!is.null(input$drug_class) && length(input$drug_class) > 0) {
      filtered <- filtered %>% filter(drug_class %in% input$drug_class)
    }
    if (!is.null(input$phase_selection) && length(input$phase_selection) > 0) {
      filtered <- filtered %>% filter(Reaction_phases %in% input$phase_selection)
    }
    if (!is.null(input$reaction_simplified) && length(input$reaction_simplified) > 0) {
      filtered <- filtered %>% filter(reaction_simplified %in% input$reaction_simplified)
    }
    
    filtered <- filtered %>%
      mutate(
        mz = sprintf("%.4f", mz),
        monoisotopic.Mass = sprintf("%.4f", monoisotopic.Mass)
      )
    return(filtered)
  })
  
  # LC-MS Based Search: Filter Data
  filtered_data_lcms <- eventReactive(input$lcms_search_button, {
    filtered <- expanded_data
    if (!is.na(input$mz_value) && !is.na(input$ppm_tolerance) && input$ppm_tolerance > 0) {
      mz_value <- as.numeric(input$mz_value)
      ppm_tolerance <- as.numeric(input$ppm_tolerance)
      mz_tolerance <- mz_value * ppm_tolerance / 1e6
      filtered <- filtered %>% filter(abs(mz - mz_value) <= mz_tolerance)
    }
    if (!is.na(input$rt_value)) {
      filtered <- filtered %>% filter(abs(rt - input$rt_value) <= 30)
    }
    if (!is.null(input$adduct_search) && length(input$adduct_search) > 0) {
      filtered <- filtered %>% filter(adduct %in% input$adduct_search)
    }
    
    filtered <- filtered %>%
      mutate(
        mz = sprintf("%.4f", mz),
        monoisotopic.Mass = sprintf("%.4f", monoisotopic.Mass)
      )
    return(filtered)
  })
  
  # Render the Drug-Based Results Table
  output$drug_table <- renderDT({
    dt_data <- filtered_data_drug()
    if (is.null(dt_data) || nrow(dt_data) == 0) {
      return(datatable(data.frame(Message = "No results found."), options = list(dom = 't')))
    }
    
    dt_data <- dt_data %>%
      mutate(
        Enzymes = ifelse(!is.na(enzymes) & enzymes != "",
                         paste0('<a href="https://www.kegg.jp/entry/', gsub("EC ", "", enzymes),
                                '" target="_blank">', gsub("EC ", "", enzymes), '</a>'),
                         "Unknown"
        ),
        SMILES_Link = ifelse(!is.na(SMILES) & SMILES != "",
                             paste0('<a href="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/',
                                    SMILES, '/PNG" target="_blank">', SMILES, '</a>'),
                             NA
        )
      ) %>%
      dplyr::select(
        precursorID,
        name,
        reaction = reaction_simplified,
        predicted_mass = monoisotopic.Mass,
        detected,
        mz,
        detected_rt = rt,
        Reaction_phases,
        detected_adduct = adduct,
        detected_mode = mode,
        molecular_formula,
        reaction_detail,
        biosystem,
        enzymes = Enzymes,
        InChI,
        InChIKey,
        SMILES = SMILES_Link,
        drug_group,
        drug_class
      )
    
    datatable(
      dt_data,
      colnames = c("PrecursorID", "Drug Name", "Reaction", "Predicted Mass", "Detected",
                   "Detected <i>m/z</i>", "Detected RT (sec)", "Reaction Phase", "Detected Adduct",
                   "Detected Mode", "Molecular Formula", "Reaction Detail", "Biosystem",
                   "Enzymes (KEGG)", "InChI", "InChIKey", "SMILES", "Drug Group", "Drug Class"),
      rownames = FALSE,
      selection = "none",
      escape = FALSE,
      extensions = c('Responsive'),
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        responsive = TRUE,
        columnDefs = list(
          list(visible = FALSE, targets = 10:18)
        ),
        rowCallback = JS("
          function(row, data) {
            var details = '<tr class=\"details-row\"><td colspan=\"10\">' +
              '<strong>Molecular Formula:</strong> ' + data[10] + '<br>' +
              '<strong>Reaction Detail:</strong> ' + data[11] + '<br>' +
              '<strong>Biosystem:</strong> ' + data[12] + '<br>' +
              '<strong>Enzymes (KEGG):</strong> ' + data[13] + '<br>' +
              '<strong>InChI:</strong> ' + data[14] + '<br>' +
              '<strong>InChIKey:</strong> ' + data[15] + '<br>' +
              '<strong>SMILES:</strong> ' + data[16] + '<br>' +
              '<strong>Drug Group:</strong> ' + data[17] + '<br>' +
              '<strong>Drug Class:</strong> ' + data[18] +
              '</td></tr>';
            $(row).on('click', function() {
              if ($(row).next().hasClass('details-row')) {
                $(row).next().remove();
              } else {
                $(row).after(details);
              }
            });
          }
        ")
      )
    )
  })
  
  # Render the LC-MS Results Table
  output$lcms_table <- renderDT({
    dt_data <- filtered_data_lcms()
    if (is.null(dt_data) || nrow(dt_data) == 0) {
      return(datatable(data.frame(Message = "No results found."), options = list(dom = 't')))
    }
    
    dt_data <- dt_data %>%
      mutate(
        Enzymes = ifelse(!is.na(enzymes) & enzymes != "",
                         paste0('<a href="https://www.kegg.jp/entry/', gsub("EC ", "", enzymes),
                                '" target="_blank">', gsub("EC ", "", enzymes), '</a>'),
                         "Unknown"
        ),
        SMILES_Link = ifelse(!is.na(SMILES) & SMILES != "",
                             paste0('<a href="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/',
                                    SMILES, '/PNG" target="_blank">', SMILES, '</a>'),
                             NA
        )
      ) %>%
      dplyr::select(
        precursorID,
        name,
        reaction = reaction_simplified,
        predicted_mass = monoisotopic.Mass,
        detected,
        mz,
        detected_rt = rt,
        detected_adduct = adduct,
        detected_mode = mode,
        molecular_formula,
        reaction_detail,
        biosystem,
        enzymes = Enzymes,
        InChI,
        InChIKey,
        SMILES = SMILES_Link,
        drug_group,
        drug_class
      )
    
    datatable(
      dt_data,
      colnames = c("PrecursorID", "Drug Name", "Reaction", "Predicted Mass", "Detected",
                   "Detected <i>m/z</i>", "Detected RT (sec)", "Detected Adduct", "Detected Mode",
                   "Molecular Formula", "Reaction Detail", "Biosystem", "Enzymes (KEGG)", "InChI",
                   "InChIKey", "SMILES", "Drug Group", "Drug Class"),
      rownames = FALSE,
      selection = "none",
      escape = FALSE,
      extensions = c('Responsive'),
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        responsive = TRUE,
        columnDefs = list(
          list(visible = FALSE, targets = 9:17)
        ),
        rowCallback = JS("
          function(row, data) {
            var details = '<tr class=\"details-row\"><td colspan=\"9\">' +
              '<strong>Molecular Formula:</strong> ' + data[9] + '<br>' +
              '<strong>Reaction Detail:</strong> ' + data[10] + '<br>' +
              '<strong>Biosystem:</strong> ' + data[11] + '<br>' +
              '<strong>Enzymes (KEGG):</strong> ' + data[12] + '<br>' +
              '<strong>InChI:</strong> ' + data[13] + '<br>' +
              '<strong>InChIKey:</strong> ' + data[14] + '<br>' +
              '<strong>SMILES:</strong> ' + data[15] + '<br>' +
              '<strong>Drug Group:</strong> ' + data[16] + '<br>' +
              '<strong>Drug Class:</strong> ' + data[17] +
              '</td></tr>';
            $(row).on('click', function() {
              if ($(row).next().hasClass('details-row')) {
                $(row).next().remove();
              } else {
                $(row).after(details);
              }
            });
          }
        ")
      )
    )
  })
}

shinyApp(ui = ui, server = server)
