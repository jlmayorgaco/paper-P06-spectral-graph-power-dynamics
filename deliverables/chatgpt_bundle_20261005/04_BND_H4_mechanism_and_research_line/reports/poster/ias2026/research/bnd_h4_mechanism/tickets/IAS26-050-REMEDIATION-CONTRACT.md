# IAS26-050 — contrato de intervención

**Pregunta:** ¿puede la acción física registrada mejorar el margen transversal completo sin seleccionarse mirando el holdout?

**Tratamiento primario congelado:** H4 all-GFL, `g=0.03625 → 0.25`; red, loads, dispatch, ubicación, rating y demás controles compartidos. Escenario sin cambio como control emparejado. Un comparador adicional sólo entra con biblioteca y selección fijadas sobre datos de desarrollo antes del holdout.

**Antes de MC:** verificar la reproducción nominal del tratamiento en el DAE completo y guardar predicción de `Delta alpha`, regla de factibilidad, umbral `tau_dec`, modo/familia, coste si se usa la palabra «mejor», y manifiesto de parámetros. Si la predicción usa raíz reducida, exige M1 PASS; la comparación DAE directa puede seguir sin M1.

**Gate:** acción y reglas de evaluación publicadas en config versionada antes de IAS26-060. Si el nominal no se reproduce, marcar intervención no validada y no prometer rescate en el póster.
