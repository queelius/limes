# LimesCompilerOptions.cmake
# Compiler options module for the limes library

function(limes_apply_compiler_options target)
    # Platform-specific optimizations
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} INTERFACE
            $<$<CONFIG:Release>:-O3 -march=native -ffast-math>
            $<$<CONFIG:Debug>:-O0 -g -fsanitize=address,undefined>
        )
        target_link_options(${target} INTERFACE
            $<$<CONFIG:Debug>:-fsanitize=address,undefined>
        )

        if(ENABLE_AVX2)
            target_compile_options(${target} INTERFACE -mavx2 -mfma)
        endif()
    endif()

    if(MSVC)
        target_compile_options(${target} INTERFACE
            $<$<CONFIG:Release>:/O2 /arch:AVX2>
            $<$<CONFIG:Debug>:/Od /MDd>
        )
    endif()
endfunction()
