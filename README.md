# lambdastar — Messkern, Modelladapter und Bootstrap

Eigenständiges R-Paket für den Mess- und Projektionskern von Lambda Star,
das numerische Parsimony Functional und Personen-/Cluster-Bootstrap mit direkten
Perzentilintervallen. Version **0.8.2** verwendet standardmäßig den
Manuskriptschätzer für Punkte und sämtliche Ziehungen. Der alternative
Populationsschätzer ist ausdrücklich wählbar. Nuisance-Effekte werden gemeinsam
berücksichtigt. Alle Intervalle verwenden `ci.method="percentile"`.
Der [Modellraum-Lernpfad](inst/doc/MODEL_SPACE_WORKFLOW.md),
weitere Anwenderanleitungen und ausführbare Beispiele werden mitinstalliert.

Methodischer Hintergrund: Mike Hammes (2026),
[*From Null Rejection to Structured Localization: Lambda Star and the Thermodynamics of Frequentist Statistics*](https://doi.org/10.5281/zenodo.22962377)
(Preprint, DOI: 10.5281/zenodo.22962377).

## Installation und automatisierte Tests

Die veröffentlichte Version ist auf CRAN verfügbar:

```r
install.packages("lambdastar")
```

Der [Quellcode auf GitHub](https://github.com/docmh/lambdastar) enthält das
veröffentlichte Paket und seine Tests. [GitHub Actions](https://github.com/docmh/lambdastar/actions)
prüft macOS auf Apple Silicon mit R release, oldrel und devel sowie Linux
mit R release und devel und Windows mit R release. Die macOS-Jobs prüfen,
dass R nativ auf ARM64 läuft. Check-Logs und Plattformdetails sind als
Artefakte verfügbar. Die CI-Prüfung verwendet `--as-cran --no-manual`;
die PDF-Handbuchprüfung erfolgt zusätzlich vor einer CRAN-Einreichung.

## Lokal installieren

Voraussetzung ist R >= 4.1 mit `stats` und dem empfohlenen R-Paket `boot`.
Ein Source-Archiv kann mit `R CMD build .` im Paketverzeichnis erzeugt werden:

```r
install.packages("lambdastar_0.8.2.tar.gz", repos = NULL, type = "source")
library(lambdastar)
packageVersion("lambdastar")
```

Die Installation verwendet die erste beschreibbare Standard-Library beziehungsweise
die von R angebotene persönliche Library. Mit `lib = "/pfad/zur/library"` lässt
sich ein anderes Ziel wählen. Das
Paket ist unter [Apache 2.0](inst/doc/LICENSE-2.0.txt) lizenziert.
Copyright 2026 Mike Hammes. Kontakt: <mike.hammes@mikehammes.name>.
Eine Kopie der Lizenz wird unter `system.file("doc", "LICENSE-2.0.txt",
package="lambdastar")` mitinstalliert.

Ein eigenständiges Beispiel und die Lernanleitung werden mitinstalliert:

```r
source(system.file("examples", "quickstart.R", package = "lambdastar"))
system.file("doc", "LEARN_LAMBDASTAR.md", package = "lambdastar")
help("lambdastar")
```

## Lambda Star Schritt für Schritt lernen

Das [ausführbare Tutorial](inst/doc/LEARN_LAMBDASTAR.md) erklärt Messrauschen,
λ*, die Lokalisation in X_H, tatsächliche Bootstrap-Perzentile und F_H anhand
eines durchgehenden Beispiels mit bekannter Struktur. Sieben Abbildungen,
Vorhersagefragen und veränderbare Parameter machen den Rechenweg nachvollziehbar.
Eine [bereits ausgeführte Fassung](inst/tutorial/WALKTHROUGH.md) mit den
Standardwerten ist ebenfalls enthalten.

```r
source(system.file("examples", "learn_lambdastar.R", package="lambdastar"))
lesson <- lambdastar_tutorial(output_dir="my-lambda-star-lesson")
```

In interaktivem R hält der Ablauf vor jeder Auswertung für eine Vorhersage an.
`pause=FALSE, show_plots=FALSE` führt ihn ohne Dialog und Grafikfenster aus.
Der Ausgabeordner enthält den bebilderten `WALKTHROUGH.md`, CSV-Tabellen,
synthetische Daten, Bootstrap-Ziehungen und `lesson.rds`. Bestehende Ergebnisse
werden nicht überschrieben. Default: B=5000 und 95%; für kurze Übungsläufe
kann beispielsweise `B=199` ausdrücklich gesetzt werden.

[Fortgeschrittene Beispiele](inst/doc/ADVANCED_EXAMPLES.md) behandeln separat
Populationsschätzung, Standortadjustierung, Schulcluster und die beiden Nullfälle.
Der [Modellraum-Lernpfad](inst/doc/MODEL_SPACE_WORKFLOW.md) führt anschließend
von der Forschungsfrage zu X_H, bedingten Effekten, Kontrasten, Basisrezepten und
F_H-Vergleichen. Seine Scope-Tabelle grenzt die derzeit unterstützten Hypothesen
und Referenzräume ab. Ein eigenes Beispiel zeigt verfügbare und undefinierte
Intervalle auf zwei synthetischen E4-Datensätzen:

```r
source(system.file("examples", "intervals_and_targets.R", package="lambdastar"))
intervals <- lambdastar_intervals_example()  # B=5000, 95%, seed=1
intervals$intervals
intervals$fits$weak$bootstrap$draw_reasons
```

## Personen über mehrere Messzeitpunkte

`lambda_panel()` ergänzt feste Personen-/Zeitpunktmodelle. Jede Zeile enthält
eine Person an einem Zeitpunkt und mindestens zwei parallele Messindikatoren.
Personenintercepts werden standardmäßig durch Gruppenzentrierung adjustiert;
weitere Ziel- und Adjustierungsterme kommen aus dem deklarierten Modell.
Der Bootstrap zieht vollständige Personenverläufe mit eigenen Identitäten für
mehrfach gezogene Kopien. [Herleitung und Lernpfad](inst/doc/PANEL_MODELS.md).

```r
source(system.file("examples", "panel_models.R", package="lambdastar"))
panel <- lambdastar_panel_example()  # zwei Analysen mit B=5000, 95%, seed=1
panel$trend$estimates
confint(panel$interaction)         # Gruppe × Zeit nach Personen- und Zeitadjustierung
panel$comparison$delta_F           # Zeittrend versus Krümmung bei gleichem Nuisance-Raum
panel$trend$panel$map               # tatsächliche Personen-/Beobachtungszuordnung
```

Der Mess-Bath muss im verbleibenden Beobachtungsraum homogen und isotrop sein.
Zeitpunkte sind keine parallelen Replikate einer unveränderten Antwort.
Ungleiche Beobachtungszahlen sind möglich; explizite Complete-Case-Auswahl
entfernt ganze Personen. T_q bleibt auf der Skala des beobachteten Itemmittels.

## Mittelwert gegen Referenz und gepaarte Differenzen

Seit 0.7.0 hält `lambda_reference()` die Abweichung von einer vorgegebenen
Referenz fest, einschließlich der Mittelwertrichtung. `lambda_paired()` bildet
`after - before` für jedes parallele Itempaar und schätzt T aus diesen Differenzen.
Herleitung, Auflösungskonvention und Beispiel: [Referenzmittelwerte](inst/doc/REFERENCE_MEANS.md).

```r
source(system.file("examples", "reference_means.R", package="lambdastar"))
means <- lambdastar_reference_example()  # B=5000, 95%, seed=1
means$one_sample$estimates              # Mittelwert gegen 50
means$paired$estimates                  # mittlere Veränderung gegen 0
means$paired_reference$estimates        # Veränderung gegen .5
means$paired$mean_departure             # vorzeichenbehaftete mittlere Abweichung
confint(means$paired)
```

Bei festem Ursprung ohne Nuisance gilt d=N, `~1` beschreibt die freie
Mittelwertrichtung und `~0` den festen Nullwert. Ein Vergleich verlangt denselben
Referenzvektor und Nuisance-Raum. T_q wird ausdrücklich auf der Mittelwert- bzw.
Differenzskala vorgegeben. Auch für Paare sind mindestens zwei parallele
Messpaare notwendig; eine einzelne Vorher-/Nachher-Spalte identifiziert keinen
separaten Mess-Bath. Gewöhnliche zentrierte Aufrufe bleiben unverändert.

## Kurzes Beispiel

```r
library(lambdastar)
set.seed(17)
d <- data.frame(x = rnorm(80), z = rnorm(80),
                group = factor(rep(c("a", "b"), 40)),
                condition = factor(rep(c("a", "b", "c", "d"), 20)))
signal <- d$x + 0.5 * d$z
d[paste0("item", 1:4)] <- replicate(4, signal + rnorm(80))
items <- paste0("item", 1:4)
d$y_mean <- rowMeans(d[items])
T_q <- 1/12  # Beispielannahme auf Mittelwertskala, keine geschätzte Auflösung

measurement <- lambda_measure(d, items)
two_group <- lambda_model(d, ~ group, items, T_q)       # Gruppenmittelwert-Hypothese
anova <- lambda_model(d, ~ condition, items, T_q)       # Gesamtmodell
ancova <- lambda_model(d, ~ condition + x, items, T_q)
regression <- lambda_model(d, y_mean ~ x + z, items, T_q)
interaction <- lambda_model(d, ~ group * x, items, T_q)
summary(regression)
regression$parsimony[c("rss", "dimension_cost", "geometry_cost", "free_energy")]

# Dieselbe Kodierung über explizites X_H oder ein geprüftes lm-Objekt:
via_matrix <- lambda_model(d, as.matrix(d[c("x", "z")]), items, T_q)
via_lm <- lambda_model(d, lm(y_mean ~ x + z, d), items, T_q)
stopifnot(all.equal(regression$parsimony$free_energy,
                    via_lm$parsimony$free_energy))

# Ein Datenobjekt, zwei Modelle, zweimal dieselbe Einzelmodellberechnung:
comparison <- lambda_compare(d, ~ x, ~ x + z, items, T_q)
print(comparison)
stopifnot(all.equal(comparison$delta_F,
  comparison$model_b$parsimony$free_energy - comparison$model_a$parsimony$free_energy))

# Derselbe Raum, andere Kodierung: Projektion gleich, Geometrie ggf. verschieden.
coding <- lambda_compare(d, ~ condition, ~ condition, items, T_q,
  contrasts_a = list(condition = "contr.treatment"),
  contrasts_b = list(condition = "contr.helmert"))
print(coding)
```

Items beziehungsweise Wiederholungen müssen bereits korrekt kodiert sein und
dieselbe stabile Struktur messen. Der Kern nimmt parallele Messungen mit
unabhängigen, homogenen Messfehlern an. Er prüft diese Messannahmen nicht
empirisch und interpretiert zeitliche Veränderungen nicht automatisch als Rauschen.
Die ausgewählten Eingaben müssen vollständig und endlich sein. Ein Ausschluss
fehlender Fälle muss ausdrücklich angefordert werden.

## Modelladapter und gemeinsame Fälle

- `lambda_measure(data, indicators)` schätzt den Messzustand. Nur hier bedeutet
  `indicators = NULL`, dass alle Spalten Messungen sind.
- `lambda_model(data, model, indicators, T_q)` wertet ein ganzes Modell aus.
  `model` ist ein `lambda_design`, eine Formel, ein geeignetes `lm` oder eine numerische Matrix.
  Einseitige Formeln verwenden den Itemmittelwert. Bei zweiseitigen Formeln
  muss die Antwort damit übereinstimmen; Abweichungen lösen einen Fehler aus.
- `lambda_compare(data, model_a, model_b, indicators, T_q)` bereitet gemeinsame
  Daten und Messbedingungen einmal vor und ruft denselben Auswertungsweg zweimal
  auf. `delta_F = F_Hb - F_Ha`; ein negatives Delta bedeutet kleineres F_H für b.
  Das ist kein Signifikanztest und keine automatische Modellauswahl.

`na.action = "fail"` ist der Default. `na.action = "complete"` wählt einmal
vollständige Personen über Items und alle Variablen beider Modelle aus.
Nicht benötigte Datenspalten beeinflussen die Auswahl nicht. Ausgeschlossene
Positionen stehen in `$cases$excluded`, verwendete in `$cases$included`;
Zeilenidentitäten und benötigte Variablen werden ebenfalls gespeichert.
Zwei getrennte Einzelaufrufe sind bei Missingness nur dann vergleichbar, wenn
beide auf genau dieser gemeinsamen Teilstichprobe ausgewertet werden.

Die Fallauswahl erfolgt **vor** Transformationen wie `poly()` oder `scale()`.
Durch Transformationen entstehende NA/Inf werden abgewiesen, nicht zusätzlich
pro Modell herausgefiltert. Formelvariablen müssen aus `data` stammen.
Eine explizite Matrix enthält alle Eingabezeilen in Datenreihenfolge; vorhandene
Zeilennamen müssen übereinstimmen. Ohne Zeilennamen wird die angegebene
Reihenfolge verwendet. Fehlende Matrixwerte gehen in die gemeinsame Auswahl ein.

`lm` benötigt einen gespeicherten Modellrahmen, einen Intercept, dieselben Fälle
in derselben Reihenfolge, eine passende Mittelwertantwort und passende kodierte
Prädiktorwerte. Gespeicherte Kontraste gelten auch nach Änderungen der globalen
R-Kontrastoptionen. Gewichte, Offsets, GLM, multivariate Antworten und andere
Modellklassen werden abgewiesen. Abweichende `lm`-Stichproben werden ausdrücklich
neu angepasst; der Adapter führt keinen stillen Refit aus.

Ein Intercept zählt nicht als Hypothesenrichtung. `~ 1` ist das Nullmodell.
Faktorstufen und tatsächliche Kontraste bleiben erhalten. Explizite
Kontrastmatrizen dürfen weniger als G−1 Spalten haben; genau diese Richtungen
bleiben auch im Bootstrap erhalten. Rangmangel, verlorene
Stufen oder ein nach Zentrierung redundanter Dummyblock wie `~ 0 + factor`
führen zum Fehler, statt die Kodierung still zu ändern. `$model` enthält Formel,
Spalten, Termzuordnung, Kontraste, Faktorstufen, Rang und Fälle; `$basis$X`
enthält das zentrierte, gegebenenfalls nuisance-adjustierte kodierte Design.
Nur die Projektion nutzt die QR-Basis.

Der Manuskriptpfad setzt stabile Struktur plus homogenen isotropen Mess-Bath
mit endlichen zweiten Momenten voraus. Weder Normalität der Antwort noch gleiche
Gesamtvarianzen der Gruppen sind erforderlich. Gruppenvergleiche (auch neben
Welch), feste ANOVA/ANCOVA, Regression und Interaktionen verwenden dasselbe X_H.
Welch-Teststatistiken werden separat berechnet; sie ersetzen weder T noch die
Projektionskorrektur. Der Populationsschätzer behält seine zusätzlichen Annahmen.

`lambda_design()` bindet ein wiederverwendbares Modell an Kodierung und
Eingabezeilen. `lambda_hypothesis()` wählt Ziel- und Adjustierungsterme aus dem
vollständigen Modell. Bedingte λ* und κ beziehen sich auf die nach Adjustierung
verbleibende Struktur. Anleitung: [Modellräume](inst/doc/MODEL_SPACES.md).

```r
spec <- lambda_design(d, ~ group * x)
effect <- lambda_hypothesis(d, ~ group * x, target="group:x", adjust=c("group", "x"))
conditional <- lambda_model(d, effect, items, T_q)
```

`lambda_contrast()` kodiert benannte homogene Restriktionen `L beta=0`;
`lambda_marginal_contrast()` verwendet ausdrücklich gesetzte Prädiktorwerte und
Gewichte für einfache oder gemittelte Effekte. `encoding="fixed"` hält die nach
Fallauswahl erlernte Basis während des Bootstrap fest; `"reevaluate"` bleibt der
Default. Anleitung und ausführbares Beispiel: [Kontraste und Basisrezepte](inst/doc/CONTRASTS_AND_BASES.md).

```r
h <- lambda_contrast(d, ~group*x, c(groupb=1, `groupb:x`=1))
fit <- lambda_model(d, h, items, T_q)  # b-a bei x=1, bedingt auf das Nullmodell
source(system.file("examples", "contrasts_and_bases.R", package="lambdastar"))
lesson <- lambdastar_contrasts_example(B=199)  # kurzer Übungslauf; Default B=5000
```

Ein-Stichproben-/gepaarte Mittelwerthypothesen verwenden die explizite
[Referenz-API](inst/doc/REFERENCE_MEANS.md). Die gewöhnliche Zentrierung entfernt
weiterhin den Mittelwert. Feste Personen-/Zeitpunktmodelle verwenden
[lambda_panel](inst/doc/PANEL_MODELS.md); allgemeine Kovarianz- und Mixed Models
bleiben gesonderte Erweiterungen.

`measurement_scale = "mean"` ist der Default; `"item"` ändert nur den
Nonzentralitätsnenner. Der Adapter verwendet für F_H in beiden Fällen korrekt
die Mittelwertantwort und `temperature_mean`. `T_q` ist dafür explizit positiv
anzugeben. `capacity`, `numerical_control`, `eta_lookup` und `response_tol` sind
überschreibbar und werden in `$settings` gespeichert. Die Antwort-/lm-Prüfung
verwendet standardmäßig punktweise
`abs(a-b) <= sqrt(.Machine$double.eps) * max(1, abs(a), abs(b))`.
`print()` und `summary()` zeigen Punktschätzungen und die F_H-Zerlegung.

## Perzentilintervalle und gemeinsame Personenziehungen

`bootstrap=TRUE` verwendet `ci.method="percentile"`, B=5000, 95 % und Seed=1.
Der Bootstrap zieht N vollständige Personen mit Zurücklegen. Jede Ziehung
berechnet T, λ*, λ_H, λ_perp und κ_H mit demselben ausgewählten Schätzer wie die
Originaldaten (`estimator="paper"` als Default). Für 95 % werden direkt die 2,5%- und 97,5%-Perzentile verwendet,
berechnet mit `boot::boot.ci(type="perc")`. Konstante Verteilungen behalten
identische untere und obere Grenzen. Niveau, B und Seed sind überschreibbar.

Die **ausgewählten Punktschätzungen stehen in `$estimates`** und in den
entsprechenden Feldern `$lambda_star`, `$lambda_h` und `$kappa`, auch bei
`bootstrap=FALSE`. Die CI-Tabelle enthält dieselben Punkte. `$inference` und
`$settings$estimator` benennen die Methode. `$paper` bewahrt die
Manuskriptwerte; `$population` wird nur bei ausdrücklicher Auswahl angelegt.

```r
ci_model <- lambda_model(d, ~ x + z, items, T_q,
  bootstrap = lambda_bootstrap_control(B = 499, conf.level = .95, seed = 17,
                                       keep_draws = TRUE))
ci_model$estimates
confint(ci_model)
confint(ci_model, parm = c("lambda_star", "lambda_h", "kappa"), level = .90)
ci_model$bootstrap$ci

# Gemeinsame Personenziehungen, feste beobachtete Gruppenanteile:
ci_comparison <- lambda_compare(d, ~ group, ~ group + x + z, items, T_q,
  bootstrap = lambda_bootstrap_control(B = 499, strata = "group"))
confint(ci_comparison)
```

`estimator="paper"` verwendet die bath-korrigierten Projektionsenergien des
Manuskripts: auf der Mittelwertskala λ_H,raw=Q_H/T−p und
λ*,raw=Q_gesamt/T−d. Ohne weitere Nuisance-Variablen ist d=N−1.
Originalschätzung und jede Bootstrap-Ziehung verwenden diese Formeln.

`estimator="population"` wählt ausdrücklich die alternative Korrektur für
Population bei festem N: λ*ₙ=(N−1)Var(S)/T und λ_H,ₙ aus der linearen
Populationsprojektion. Sie verwendet `V_S=Var(Y_mean)-T_mean` und
`V_H=Var(Y_mean)-RSS/(N-1-p)` unter linearer bedingter Erwartung und homogener
Residualvarianz. Sie unterstützt derzeit Intercept-Adjustierung und
Personen-Bootstrap. Nuisance- oder Cluster-Kombinationen werden bei dieser
Option ausdrücklich abgelehnt. `help("lambda_inference")` enthält die Formeln.

### Nuisance-Effekte und Schulcluster wie im Manuskript

```r
# sses: eine Zeile je Person, fünf bereits aufbereitete WHO-5-Items.
# skills: Namen der 15 Prädiktoren; Originalauflösung des Itemmittels.
fit_sses <- lambda_model(sses, as.matrix(sses[skills]), items, T_q = 1/300,
  nuisance = ~ factor(SiteID),
  bootstrap = lambda_bootstrap_control(B = 1000, seed = 2026073116,
    resampling = "cluster", cluster = "school_uid", strata = "SiteID",
    percentile_interpolation = "type7", keep_draws = TRUE))
fit_sses$estimates
confint(fit_sses, parm = "kappa")
```

`nuisance` akzeptiert eine einseitige Formel oder eine numerische Matrix.
Ein Intercept gehört stets zum Nuisance-Raum. Messungen und Hypothesenmatrix
werden gemeinsam residualisiert; d=N−rank(B) bestimmt die Messkorrektur.
Der Nuisance-Raum wird in jeder Ziehung neu geschätzt. Kodierung und Rangverlust
werden geprüft. Die Modellformel enthält die zu prüfenden Richtungen,
Nuisance-Variablen stehen separat in `nuisance`.

Cluster werden vollständig und mit Zurücklegen innerhalb ihrer Strata gezogen;
je Stratum bleibt die Clusterzahl fest, die Personenzahl kann variieren.
Clusterbezeichnungen gelten innerhalb des jeweiligen Stratums. Die gemeinsame
Fallauswahl umfasst auch Nuisance-, Cluster- und Stratumvariablen.
`percentile_interpolation="type7"` berechnet die unveränderten R-Perzentile
(`quantile(type=7)`) der Originalanalyse. Default `"boot"` nutzt weiterhin
`boot::boot.ci(type="perc")`; beide sind gewöhnliche Perzentilintervalle.

| Situation | Behandlung |
|---|---|
| λ*=0 | gültiger Nullwert; bleibt in der λ*-Verteilung |
| κ bei geschätztem Gesamtsignal null | `NA`, Grund `undefined_kappa` |
| κ=0 oder κ=1 bei positivem Gesamtsignal | gültige Werte; bleiben in der Verteilung |
| Originaldaten mit T=0 | Singularität; keine Intervalle |
| Bootstrap-Ziehung mit T=0 | T=0 bleibt erhalten; davon abhängige Kennwerte sind undefiniert |
| Rangverlust oder fehlende Faktorstufe | betroffene Modellwerte `NA`, eigener Diagnosegrund |

Ein undefinierter Originalschätzer erhält `undefined_estimate`. Sobald eine
Ziehung für einen Kennwert undefiniert ist, bleibt dessen reguläres Intervall
`NA` (`undefined_draws` bei κ ohne Signal, sonst `failed_draws`). Die Intervalle
anderer Kennwerte werden unabhängig davon berechnet. `n_valid`, `n_failed`,
`n_undefined` und die Ausfallgründe sind sichtbar; Ziehungen werden nicht ersetzt.
Mindestens `min_valid=200` Werte werden für Quantile benötigt.

Mit `conditional_quantiles=TRUE` können zusätzlich Quantile ausschließlich der
gültigen Ziehungen angefordert werden. Sie stehen separat in
`$bootstrap$conditional_quantiles`, mit Zahl der ausgeschlossenen Ziehungen,
als **bedingte beschreibende Quantile**. `confint()` liefert weiterhin nur die
regulären Konfidenzintervalle.

`strata="group"` hält beim Personen-Bootstrap die beobachteten Gruppengrößen fest.
Nur `estimator="population"` verwendet dabei die korrigierte Mischungsvarianz
und benötigt mindestens zwei Personen je Stratum. Ein gesättigtes Modell hat
keine Residualfreiheitsgrade für die Populationsprojektion; seine betroffenen
Populationswerte werden mit Begründung als `NA` ausgegeben.

`keep_draws=TRUE` speichert **die tatsächlichen Kennwerte** jeder Ziehung,
einschließlich null und NA, sowie `draw_reasons`. `keep_indices=TRUE` speichert
Originalzeilenpositionen, bei Clustern als Liste variabler Länge. `confint(..., level=...)` verwendet vorhandene Ziehungen
und benötigt dafür `keep_draws=TRUE`. Alle Kennwerte und Vergleichsmodelle teilen
dieselben Personenziehungen. Die Intervalle sind punktweise; F_H und ΔF bleiben
Punktschätzungen. Die Validierung berichtet Verfügbarkeit, Breite und Überdeckung
separat; bei wahrem Nullsignal ist eine Überdeckung für κ nicht definiert.

## Funktionen und Skalen

- `estimate_lambda_star_items()` liefert Temperatur, rohe und nichtnegative
  Gesamtstruktur sowie Zeilenidentität. Default ist `scale = "composite_mean"`.
- `estimate_lambda_partition_items()` ergänzt die Projektion auf ein explizites
  `x`, die rohe/beschränkte Zerlegung und Kappa. Das Design enthält keinen Intercept.
- `hypothesis_basis()` bewahrt die kodierte Matrix; die normalisierte QR-Basis
  dient nur der Projektion. `design_volume(logarithm = TRUE)` bewahrt log(V)
  auch bei einem sehr kleinen Volumen.
- `eta_numerical()` löst die quantisierte Kapazität über eine logarithmische
  Momentengleichung. `details = TRUE` zeigt Lösungs- und Abschneidediagnostik.
- `xi_eta_numerical()` integriert standardmäßig direkt und adaptiv. Eine
  wiederverwendbare Tabelle aus `make_eta_numerical_lookup()` beschleunigt viele
  Auswertungen; unterhalb ihres ersten positiven Knotens wird direkt integriert.
- `thermodynamic_model_state()` berechnet F_H für ein einzelnes X_H. Der
  Default ist `eta_method = "numerical"`; `"smooth"` ist eine ausdrücklich
  gekennzeichnete Approximation. Beide verwenden dieselbe geometrische Aktivierung.
- `lambda_control()` hält überschreibbare Toleranzen getrennt für Rang,
  Kapazitätslösung, abgeschnittene Verteilungsschwänze, Quadratur und Lookup.
  Lookup-Fehler sind geschätzte Diskretisierungsfehler, keine formalen Schranken
  für den gesamten Rechenweg. Die Gibbs-Tails werden durch Integralgrenzen geprüft.

`scale = "item"` ändert den Nenner der Nonzentralitäten auf die
Einzelmessungsvarianz. Das gespeicherte `composite_response` und die
`stable_profile` sind weiterhin Mittelwerte. Wer die Mittelwertantwort in F_H
verwendet, muss deshalb **temperature_mean** verwenden. Eine Antwortskalierung
um den Faktor c erfordert die entsprechende Skalierung von T und T_q um c².
Die Auflösung T_q wird fachlich vorgegeben.

Bei T = 0 liefert der Kern `status = "singular_temperature"` und NA für die
nicht regulär schätzbaren normierten Größen beziehungsweise F_H. Für explizite
numerische Grenzwertprüfungen ist `boundary = TRUE` möglich: F_H = RSS,
Status `boundary_only`. Dies ist kein empirischer Schätz- oder Vergleichszustand.
Kappa bleibt auch bei positiver Temperatur und geschätzter Gesamtstruktur null
nicht bestimmbar. Sehr kleine positive Werte werden nicht durch epsilon ersetzt.

## Bauen und prüfen

Das Paket importiert `stats` und `boot`. Der Testrunner verwendet Base R;
Dokumentation und NAMESPACE sind bereits enthalten.

```sh
R CMD build .
R CMD check --as-cran lambdastar_0.8.2.tar.gz
```

Die Tests enthalten eine vorhandene synthetische Simulation und unabhängig
mit Python Decimal erzeugte numerische Referenzwerte. Die Anwenderbeispiele
verwenden synthetische Daten; externe OECD-Rohdaten werden nicht mitgeliefert.
Die Apache-2.0-Lizenz gilt für dieses Paket; sie lizenziert weder das separate
Manuskript noch externe OECD-Daten um.
