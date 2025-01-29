
/*
Este codigo de testing es ideal para probar el archivo de configuración y además hacer pruebas si existen probelmas de escreibir datos
permitiendo mayor cantidad de debugging
*/
int main()
{
    try
    {
        // 1. Leer archivo de configuración
        Config cfg = readConfig("parametros.conf");

        // 2. Calcular total de iteraciones
        double total_iterations = computeTotalIterations(cfg);
        cout << "Total de iteraciones: " << total_iterations << endl;

        double current_iteration = 0.0;

        // 3. Recorrer los bucles usando los valores del cfg
        for(double lambda1 = cfg.lambda1_min; lambda1 <= cfg.lambda1_max; lambda1 += cfg.step_lambda1)
        {
            for(double lambda2 = cfg.lambda2_min; lambda2 <= cfg.lambda2_max; lambda2 += cfg.step_lambda2)
            {
                for(double lambda3 = cfg.lambda3_min; lambda3 <= cfg.lambda3_max; lambda3 += cfg.step_lambda3)
                {
                    for(double lambda4 = cfg.lambda4_min; lambda4 <= cfg.lambda4_max; lambda4 += cfg.step_lambda4)
                    {
                        for(double lambda5 = cfg.lambda5_min; lambda5 <= cfg.lambda5_max; lambda5 += cfg.step_lambda5)
                        {
                            for(double m12_squared = cfg.m12_squared_min; 
                                        m12_squared <= cfg.m12_squared_max; 
                                        m12_squared += cfg.step_m12_squared)
                            {
                                for(double beta = cfg.beta_min; beta <= cfg.beta_max; beta += cfg.step_beta)
                                {
                                    // Aquí hacemos lo que necesitemos con los parámetros (lambda1, lambda2, etc.)
                                    // ...

                                    // Actualizamos el contador
                                    current_iteration++;
                                    // (Opcional) imprimir avance cada cierto número de iteraciones
                                    if (fmod(current_iteration, 100000) == 0) {
                                        cout << "Progreso: " << (current_iteration / total_iterations) * 100.0
                                             << "%\r" << flush;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        cout << "\nFinalizado. Iteraciones totales procesadas: " << current_iteration << endl;
    }
    catch(const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
