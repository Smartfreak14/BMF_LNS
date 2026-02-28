#!/bin/bash
echo "=== Build RBAC Maintenance C++ (avec EvalMaxSAT2022 / IPAMIR) ==="

# Se positionner dans le répertoire du script
cd "$(dirname "$0")"

# Créer répertoires
mkdir -p build

# Compiler EvalMaxSAT2022 si la bibliothèque n'existe pas
EVALMAXSAT2022_DIR="third_party/EvalMaxSAT2022"
EVALMAXSAT2022_LIB="${EVALMAXSAT2022_DIR}/libipamirEvalMaxSAT2022.a"

if [ ! -f "$EVALMAXSAT2022_LIB" ]; then
    if [ -d "$EVALMAXSAT2022_DIR" ]; then
        echo "Compilation d'EvalMaxSAT2022 (IPAMIR)..."
        cd "$EVALMAXSAT2022_DIR"
        
        # Compiler via makefile
        make clean 2>/dev/null
        make
        
        if [ -f "libipamirEvalMaxSAT2022.a" ]; then
            echo "✅ EvalMaxSAT2022 compilé avec succès"
        else
            echo "❌ Échec de la compilation d'EvalMaxSAT2022"
            exit 1
        fi
        
        cd ../..
    else
        echo "❌ EvalMaxSAT2022 non trouvé dans $EVALMAXSAT2022_DIR"
        echo "   Veuillez ajouter EvalMaxSAT2022 dans third_party/"
        exit 1
    fi
else
    echo "✅ EvalMaxSAT2022 déjà compilé: $EVALMAXSAT2022_LIB"
fi

# Build project principal
echo "Building main project..."
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(sysctl -n hw.ncpu 2>/dev/null || nproc)
cd ..

echo "=== Build terminé ==="
if [ -f "build/rbac_main" ]; then
    echo "✅ Succès! Usage: ./build/rbac_main"
    echo ""
    echo "Exécutables disponibles:"
    echo "  ./build/rbac_main      - Programme principal"
    echo "  ./build/test_nolk      - Tests NoLk (confidentialité)"
    echo "  ./build/test_nocorrupt - Tests NoCorrupt (intégrité)"
    echo "  ./build/test_all       - Tests complets"
else
    echo "❌ Échec du build"
fi
