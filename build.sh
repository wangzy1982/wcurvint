PROJECT_ROOT=$(cd "$(dirname "$0")" && pwd)
BUILD_DIR="${PROJECT_ROOT}/build"
BIN_DIR="${BUILD_DIR}/bin"
WSOLVER_DLL_SRC="${PROJECT_ROOT}/thirdpart/wsolver/dll/wsolver.dll"

echo "===================== 1. Create Build Directory ====================="
if [ ! -d "${BUILD_DIR}" ]; then
    mkdir -p "${BUILD_DIR}"
    echo "Successfully created build directory: ${BUILD_DIR}"
else
    echo "Build directory already exists: ${BUILD_DIR}"
fi

echo -e "\n===================== 2. CMake Configuration (Generate Build Files) ====================="
cd "${BUILD_DIR}" || {
    echo "Error: Failed to enter build directory ${BUILD_DIR}"
    exit 1
}

cmake .. -A x64
if [ $? -ne 0 ]; then
    echo "Error: CMake configuration failed. Please check CMakeLists.txt or dependent files."
    exit 1
fi
echo "CMake configuration succeeded. VS solution generated in ${BUILD_DIR}"

echo -e "\n===================== 3. Compile Project (Release Version) ====================="
cmake --build . --config Release --parallel

if [ $? -ne 0 ]; then
    echo "Error: Project compilation failed. Please check code or dependent libraries."
    exit 1
fi
echo "Project compilation succeeded. Products output to ${BIN_DIR}"

echo -e "\n===================== 4. Copy wsolver.dll to Output Directory ====================="
if [ ! -f "${WSOLVER_DLL_SRC}" ]; then
    echo "Warning: wsolver.dll not found at path: ${WSOLVER_DLL_SRC}. Skipping copy step."
    exit 0
fi

cp -fv "${WSOLVER_DLL_SRC}" "${BIN_DIR}/"

if [ $? -eq 0 ]; then
    echo "wsolver.dll copied successfully. Target path: ${BIN_DIR}/wsolver.dll"
else
    echo "Error: Failed to copy wsolver.dll."
    exit 1
fi

echo -e "\n===================== All Operations Completed Successfully! ====================="
echo "Build output directory: ${BIN_DIR}"
echo "Included files: wcurvint.dll (core library), wsolver.dll (dependent library), wcurvint.lib (import library, located in ${BUILD_DIR}/lib)"