# Volumetric Clouds

This is a demo implementing one way of drawing volumetric
cloudscapes in Godot sky shaders.

This demo uses animated clouds generated from ray marching
3D textures. It features automatic time of day adjustments
just by rotating the sun.

This is version 2. This approach is more in line with modern
approaches to drawing clouds. It also relies on some Godot
features only present in version 4.2 (i.e. TextureRD for
writing to a texture using compute shaders and reading that
texture in the Godot scene).

This demo requires Godot 4.2 or later.

2 optimizations:
- Render the hemisphere to a texture over 64 frames
    - interpolate between two copies of that texture to hide changes
- Use heirarchical step lengths to reduce number of samples while raymarching

Renderer: Vulkan

## Screenshots


