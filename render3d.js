function initRender3D( sim, container ) {
  var container;

  var renderer = null, sx=0, sy=0, sz=0, rmax = 1000;
  var mesh, zmesh, lightMesh, geometry;

  var target = { x:0, y:0, z:0 };
  var intHandler = null;
  var dragging = false;

  var cw = $(container).width();
  var ch = $(container).height();

  var camera = new THREE.Camera( 30, cw/ch, 1, 10000 ); // aspect ratio 800/500!
  var scene = new THREE.Scene();
  scene.fog = new THREE.Fog( 0xffffff, 1, 10000 );

  var geometry;

  // select renderer
  if (typeof Float32Array != "undefined") {
    try {
      renderer = new THREE.WebGLRenderer( { clearColor: 0xffffff ,clearAlpha: 1 }  );
    } catch(e) { }
  }
  if( !renderer ) {
    renderer = new THREE.CanvasRenderer( { clearColor: 0xffffff ,clearAlpha: 0 }  );
  }
  renderer.setSize( cw, ch );
  geometry = new THREE.Geometry();

  // atoms
  for( i = 0; i < sim.atoms.length; ++i )
  {
    // TODO: element styling
    vector = new THREE.Vector3( sim.atoms[i].p.x-sim.ss.x/2, sim.atoms[i].p.y-sim.ss.y/2, sim.atoms[i].p.z-sim.ss.z/2 );
    geometry.vertices.push( new THREE.Vertex( vector ) );
  }
  camera.position.x = sim.ss.x/2;
  camera.position.y = sim.ss.y/2;
  camera.position.z = -70;
  //setTarget(1000,1000,-1000)
  //geometry.colors = colors;

  var sprite = THREE.ImageUtils.loadTexture( "ball.png" );
  var material = new THREE.ParticleBasicMaterial( { size: 1.5, map: sprite, vertexColors: false } );
  material.color.setRGB( 1.0, 0.0, 0.0 );

  particles = new THREE.ParticleSystem( geometry, material );
  particles.sortParticles = true;
  particles.updateMatrix();
  scene.addObject( particles );

  var light = new THREE.DirectionalLight( 0xffffff );
  light.position.x = 1;
  light.position.y = 0;
  light.position.z = 0;
  scene.addLight( light );

  $(container).empty().append( renderer.domElement );
  /*$(renderer.domElement).bind( 'mousemove', onMouseMove );
  $(renderer.domElement).bind( 'mousedown', function() { dragging = true; } );
  $(renderer.domElement).bind( 'mouseup', function() { dragging = false; } );
  $(renderer.domElement).bind( 'mouseout', function() { dragging = false; } );*/


  function setTarget(x,y,z) {
    target.x=x;
    target.y=y;
    target.z=z;
    if( intHandler === null ) {
      intHandler = setInterval( loop, 1000 / 25 );
    }
  }

  function onMouseMove(event) {
    var radius = 2*event.clientY/ch * rmax;
    var theta  = ( event.clientX/cw - 0.5 ) * 4.0 * Math.PI;
    if( dragging ) {
      target.x = Math.sin(theta)*radius;
      target.y = 0;
      target.z = Math.cos(theta)*radius;
      if( intHandler === null ) {
        intHandler = setInterval( loop, 1000 / 25 );
      }
    }
  }

  function loop() {
    var dx = ( target.x - camera.position.x ) * .05;
    var dy = ( -target.y - camera.position.y ) * .05;
    var dz = ( target.z - camera.position.z ) * .05;

    camera.position.x += dx;
    camera.position.y += dy;
    camera.position.z += dz;
    //$('#chancoord').text( camera.position.x.toFixed() + ',' + (-camera.position.y).toFixed() );
    renderer.render( scene, camera );

    if( Math.abs(dx) < 0.1 && Math.abs(dy) < 0.1 && Math.abs(dz) < 0.1  ) {
      clearInterval(intHandler);
      intHandler = null;
    }
  }

  function updateScene() 
  {
    var i;
    // update vertices
    for( i = 0; i < sim.atoms.length; ++i )
    {
      particles.geometry.vertices[i].position.x = sim.atoms[i].p.x;
      particles.geometry.vertices[i].position.y = sim.atoms[i].p.y;
      particles.geometry.vertices[i].position.z = sim.atoms[i].p.z;
    }
    particles.geometry.__dirtyVertices = true;
    renderer.render( scene, camera );
  }

  function generateSprite() {
    var canvas = document.createElement( 'canvas' );
    canvas.width = 16;
    canvas.height = 16;

    var context = canvas.getContext( '2d' );
    var gradient = context.createRadialGradient( canvas.width / 2, canvas.height / 2, 0, canvas.width / 2, canvas.height / 2, canvas.width / 2 );
    gradient.addColorStop( 0, 'rgba(255,255,255,1)' );
    gradient.addColorStop( 0.2, 'rgba(0,255,255,1)' );
    gradient.addColorStop( 0.4, 'rgba(0,0,64,1)' );
    gradient.addColorStop( 1, 'rgba(0,0,0,1)' );

    context.fillStyle = gradient;
    context.fillRect( 0, 0, canvas.width, canvas.height );

    return canvas;
  }

  return updateScene;
};
