# RBAC Maintenance C++

Implementation C++ des algorithmes de maintenance RBAC avec préservation des propriétés de sécurité, utilisant le solveur MaxSAT **EvalMaxSAT**.

## 📋 Description

Ce projet implémente en C++ les algorithmes de maintenance RBAC (Role-Based Access Control) avec deux propriétés de sécurité:

### NoLk (Confidentialité)
- **Objectif**: Empêcher les fuites d'information entre deux sujets
- **Contrainte**: `¬(F^w[i1,j] ∧ F^r[i2,j])` pour tout objet j
- **Signification**: Le sujet i1 ne peut pas écrire sur un objet que le sujet i2 peut lire

### NoCorrupt (Intégrité)
- **Objectif**: Empêcher la corruption d'information d'un objet vers un autre
- **Contrainte**: `¬(F^r[i,j1] ∧ F^w[i,j2])` pour tout sujet i
- **Signification**: Si un sujet lit j1, il ne peut pas (même indirectement) écrire sur j2

## 🔧 Opérations Supportées

| Opération | NoLk | NoCorrupt |
|-----------|------|-----------|
| Vérification | ✅ | ✅ (SAT + Graph) |
| Add-Subject | ✅ | ✅ |
| Add-Object | ✅ | ✅ |
| Modify-Rule | ❌ | ✅ |

## 🚀 Installation

### Prérequis
- CMake >= 3.16
- Compilateur C++ avec support C++17
- Git (pour cloner EvalMaxSAT)

### Compilation

```bash
cd /Users/fotsofranck/Projets/BMF/rbac_maintenance/cpp

# Rendre le script exécutable
chmod +x build.sh

# Compiler (installe automatiquement EvalMaxSAT)
./build.sh
```

## 📖 Utilisation

### Exécutables disponibles

```bash
# Programme principal (tous les tests de base)
./build/rbac_main

# Options du programme principal
./build/rbac_main --nolk        # Tests NoLk seulement
./build/rbac_main --nocorrupt   # Tests NoCorrupt seulement
./build/rbac_main --performance # Tests de performance
./build/rbac_main --all         # Tous les tests
./build/rbac_main --help        # Aide

# Tests spécifiques
./build/test_nolk      # Tests détaillés NoLk
./build/test_nocorrupt # Tests détaillés NoCorrupt
./build/test_all       # Tests complets avec scalabilité
```

### Exemple d'utilisation en C++

```cpp
#include "RBACMaintenanceNoLk.hpp"
#include "RBACMaintenanceNoCorrupt.hpp"

int main() {
    // Définir les matrices A et B
    Matrix A({
        {1, 0},  // u0: rôle 0
        {0, 1}   // u1: rôle 1
    });
    
    Matrix B({
        {1, 0, 0, 0},  // rôle 0: read r0
        {0, 0, 1, 0}   // rôle 1: read r1
    });
    
    // ===== NoLk =====
    RBAC_Maintenance_NoLk solver_nolk(2);  // k=2 rôles
    
    // Vérifier NoLk
    auto [sat, time] = solver_nolk.verify_nolk_via_sat(A, B, 0, 1);
    
    // Ajouter un sujet
    std::vector<int> permissions = {1, 0, 1, 0};  // read r0, read r1
    auto result = solver_nolk.add_subject_nolk(permissions, 0, 1, A, B);
    
    // ===== NoCorrupt =====
    RBAC_Maintenance_NoCorrupt solver_nocor(2);
    
    // Vérifier via graphe (plus rapide)
    auto [ok, t] = solver_nocor.verify_nocorrupt_via_graph(A, B, 0, 1);
    
    // Modifier une règle
    auto mod = solver_nocor.modify_rule_nocorrupt(0, 0, "write", 0, 1, A, B);
    
    return 0;
}
```

## 📁 Structure du projet

```
cpp/
├── CMakeLists.txt          # Configuration CMake
├── build.sh                # Script de compilation
├── README.md               # Ce fichier
├── src/
│   ├── Matrix.hpp          # Classe Matrix
│   ├── SATSolver.hpp       # Interface solveur SAT/MaxSAT
│   ├── SATSolver.cpp       # Implémentation avec EvalMaxSAT
│   ├── RBACMaintenanceNoLk.hpp    # Header NoLk
│   ├── RBACMaintenanceNoLk.cpp    # Implémentation NoLk
│   ├── RBACMaintenanceNoCorrupt.hpp   # Header NoCorrupt
│   ├── RBACMaintenanceNoCorrupt.cpp   # Implémentation NoCorrupt
│   └── main.cpp            # Programme principal
├── test_nolk.cpp           # Tests NoLk
├── test_nocorrupt.cpp      # Tests NoCorrupt
├── test_all.cpp            # Tests complets
├── build/                  # Répertoire de build (généré)
└── third_party/
    └── EvalMaxSAT/         # Solveur MaxSAT (cloné automatiquement)
```

## 🔬 Correspondance Python/Cython → C++

| Python/Cython | C++ |
|--------------|-----|
| `RBAC_Maintenance_NoLk` | `RBAC_Maintenance_NoLk` |
| `RBAC_Maintenance_nocorrupt` | `RBAC_Maintenance_NoCorrupt` |
| `VariableManager` | `VariableManager` |
| `pysat.solvers.Solver` | `SATSolver` (EvalMaxSAT) |
| `pysat.examples.rc2.RC2` | `SATSolver::add_soft_clause()` |
| `flow_direct()` | `flow_direct()` |
| `flow_indirect()` | `flow_indirect()` |
| `flow_indirect_optimized()` | `flow_indirect_optimized()` |
| `nolk_constraint()` | `nolk_constraint()` |
| `nocorrupt_constraint()` | `nocorrupt_constraint()` |
| `_verify_nolk_via_sat()` | `verify_nolk_via_sat()` |
| `_verify_nocorrupt_via_graph()` | `verify_nocorrupt_via_graph()` |
| `add_subject_nolk()` | `add_subject_nolk()` |
| `add_object_nolk()` | `add_object_nolk()` |
| `add_subject_nocorrupt()` | `add_subject_nocorrupt()` |
| `add_object_nocorrupt()` | `add_object_nocorrupt()` |
| `modify_rule_nocorrupt()` | `modify_rule_nocorrupt()` |

## ⚡ Performance

Le solveur EvalMaxSAT avec CaDiCaL offre d'excellentes performances pour les problèmes MaxSAT pondérés. Les temps de résolution typiques:

| m × n_obj | k | Vérification | Ajout Sujet |
|-----------|---|--------------|-------------|
| 10 × 5 | 3 | < 1 ms | < 5 ms |
| 20 × 10 | 5 | < 5 ms | < 20 ms |
| 50 × 25 | 8 | < 50 ms | < 200 ms |
| 100 × 50 | 10 | < 200 ms | < 1 s |

## 📚 Références

- EvalMaxSAT: https://github.com/FlorentAvellaneda/EvalMaxSAT
- CaDiCaL: https://github.com/arminbiere/cadical
